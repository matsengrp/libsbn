// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "subsplit_dag.hpp"

#include "combinatorics.hpp"
#include "numerical_utils.hpp"
#include "sbn_probability.hpp"

// ** Constructor methods:

SubsplitDAG::SubsplitDAG()
    : taxon_count_(0), edge_count_without_leaf_subsplits_(0), topology_count_(0.) {}

SubsplitDAG::SubsplitDAG(const RootedTreeCollection &tree_collection)
    : SubsplitDAG(tree_collection.TaxonCount(), tree_collection.TopologyCounter(),
                  tree_collection.TagTaxonMap()) {}

SubsplitDAG::SubsplitDAG(size_t taxon_count,
                         const Node::TopologyCounter &topology_counter,
                         const TagStringMap &tag_taxon_map)
    : dag_taxa_(), taxon_count_(taxon_count) {
  Assert(topology_counter.size() > 0, "Empty topology counter given to SubsplitDAG.");
  Assert(topology_counter.begin()->first->LeafCount() == taxon_count,
         "Taxon count mismatch in SubsplitDAG constructor.");
  auto [edge_indexer, index_to_child, rootsplits] =
      ProcessTopologyCounter(topology_counter);
  BuildTaxonMap(tag_taxon_map);
  BuildNodes(index_to_child, rootsplits);
  BuildEdges(index_to_child);
  BuildDAGEdgesFromEdgeIndexer(edge_indexer);
  AddLeafSubsplitsToDAGEdgesAndParentToRange();
  StoreEdgeIds();
  CountTopologies();

  Assert(IsValidTaxonMap(),
         "SubsplitDAG::SubsplitDAG(): Given Taxon Map does not represent a valid map. "
         "Taxon ids do not cover all indices from 0 to taxon_count.");
}

SubsplitDAG::SubsplitDAG(SubsplitDAG &host_dag, HostDispatchTag)
    : storage_{host_dag.storage_, HostDispatchTag()},
      dag_taxa_{host_dag.dag_taxa_},
      subsplit_to_id_{host_dag.subsplit_to_id_},
      parent_to_child_range_{host_dag.parent_to_child_range_},
      taxon_count_{host_dag.taxon_count_},
      edge_count_without_leaf_subsplits_{host_dag.edge_count_without_leaf_subsplits_},
      topology_count_{host_dag.topology_count_},
      topology_count_below_{host_dag.topology_count_below_} {}

void SubsplitDAG::ResetHostDAG(SubsplitDAG &host_dag) {
  storage_.ResetHost(host_dag.storage_);
  dag_taxa_ = host_dag.dag_taxa_;
  subsplit_to_id_ = host_dag.subsplit_to_id_;
  parent_to_child_range_ = host_dag.parent_to_child_range_;
  taxon_count_ = host_dag.taxon_count_;
  edge_count_without_leaf_subsplits_ = host_dag.edge_count_without_leaf_subsplits_;
  topology_count_ = host_dag.topology_count_;
  topology_count_below_ = host_dag.topology_count_below_;
}

// ** Comparator

int SubsplitDAG::Compare(const SubsplitDAG &other) {
  return SubsplitDAG::Compare(*this, other);
}

int SubsplitDAG::Compare(const SubsplitDAG &lhs, const SubsplitDAG &rhs) {
  // (1) Compare Taxon Sizes.
  int taxon_diff = lhs.TaxonCount() - rhs.TaxonCount();
  if (taxon_diff != 0) {
    // #350 let's talk about -100 vs -200 here.
    return taxon_diff;
  }
  // Create translation map (lhs->rhs) for bitset clades.
  auto taxon_map = SubsplitDAG::BuildTaxonTranslationMap(lhs, rhs);
  // (2) Compare Subsplit Nodes.
  auto lhs_nodes = lhs.GetSortedVectorOfNodeBitsets();
  auto rhs_nodes = rhs.GetSortedVectorOfNodeBitsets();
  // Translate to account for different Taxon mappings and sort output.
  for (size_t i = 0; i < lhs_nodes.size(); i++) {
    lhs_nodes[i] =
        SubsplitDAG::BitsetTranslateViaTaxonTranslationMap(lhs_nodes[i], taxon_map);
    lhs_nodes[i] = lhs_nodes[i].SubsplitSortClades();
  }
  std::sort(lhs_nodes.begin(), lhs_nodes.end());
  if (lhs_nodes != rhs_nodes) {
    return (lhs_nodes < rhs_nodes) ? -1 : 1;
  }
  // (3) Compare PCSP Edges.
  auto lhs_edges = lhs.GetSortedVectorOfEdgeBitsets();
  auto rhs_edges = rhs.GetSortedVectorOfEdgeBitsets();
  // Translate to account for different Taxon mappings and sort output.
  for (size_t i = 0; i < lhs_edges.size(); i++) {
    lhs_edges[i] =
        SubsplitDAG::BitsetTranslateViaTaxonTranslationMap(lhs_edges[i], taxon_map);
    lhs_edges[i] = lhs_edges[i].PCSPSortClades();
  }
  std::sort(lhs_edges.begin(), lhs_edges.end());
  if (lhs_edges != rhs_edges) {
    return (lhs_edges < rhs_edges) ? -1 : 1;
  }
  return 0;
}

bool operator==(const SubsplitDAG &lhs, const SubsplitDAG &rhs) {
  return (SubsplitDAG::Compare(lhs, rhs) == 0);
}

bool operator!=(const SubsplitDAG &lhs, const SubsplitDAG &rhs) {
  return (SubsplitDAG::Compare(lhs, rhs) != 0);
}

// ** Count methods:

void SubsplitDAG::CountTopologies() {
  topology_count_below_ = EigenVectorXd::Ones(NodeCount());

  for (const auto &node_id : RootwardNodeTraversalTrace(true)) {
    const auto &node = GetDAGNode(node_id);
    for (const bool rotated : {true, false}) {
      // When there are no leafward nodes in the `rotated` direction, we set the number
      // of topologies for the rotation of the node to be 1.
      double per_rotated_count = node.GetLeafward(rotated).empty() ? 1. : 0.;
      // Sum options across the possible children.
      for (const auto &child_id : node.GetLeafward(rotated)) {
        per_rotated_count += topology_count_below_[child_id.value_];
      }
      // Take the product across the number of options for the left and right branches
      // of the tree.
      topology_count_below_[node_id.value_] *= per_rotated_count;
    }
  }
  topology_count_ = topology_count_below_[GetDAGRootNodeId().value_];
}

size_t SubsplitDAG::TaxonCount() const { return dag_taxa_.size(); }

size_t SubsplitDAG::NodeCount() const { return storage_.GetVertices().size(); }

size_t SubsplitDAG::NodeCountWithoutDAGRoot() const { return NodeCount() - 1; }

NodeIdPair SubsplitDAG::NodeIdRange() const { return {NodeId(0), NodeId(NodeCount())}; }

size_t SubsplitDAG::RootsplitCount() const { return GetRootsplitNodeIds().size(); }

size_t SubsplitDAG::EdgeCount() const { return edge_count_without_leaf_subsplits_; }

size_t SubsplitDAG::EdgeCountWithLeafSubsplits() const {
  return storage_.GetLines().size();
}

EdgeIdPair SubsplitDAG::EdgeIdxRange() const {
  return {EdgeId(0), EdgeId(EdgeCountWithLeafSubsplits())};
}

double SubsplitDAG::TopologyCount() const { return topology_count_; }

// ** Output methods:

StringSizeMap SubsplitDAG::SummaryStatistics() const {
  return {{"node_count", NodeCount()}, {"edge_count", EdgeCountWithLeafSubsplits()}};
}

void SubsplitDAG::Print() const {
  for (auto dag_node : storage_.GetVertices()) {
    std::cout << dag_node.ToString() << std::endl;
  }
}

void SubsplitDAG::PrintNodes() const {
  for (const auto &dag_node : storage_.GetVertices()) {
    std::cout << dag_node.Id() << ": " << dag_node.GetBitset().SubsplitToString()
              << std::endl;
  }
}

void SubsplitDAG::PrintEdgeIndexer() const {
  for (const auto &[edge, idx] : BuildEdgeIndexer()) {
    std::cout << idx << ": " << edge.PCSPToString() << std::endl;
  }
}

void SubsplitDAG::PrintDAGEdges() const {
  for (const auto &[parent_child_id, edge_idx] : storage_.GetLines()) {
    const auto &[parent_id, child_id] = parent_child_id;
    std::cout << "[" << parent_id << ", " << child_id << "], " << edge_idx << std::endl;
  }
}

void SubsplitDAG::PrintParentToRange() const {
  for (const auto &[subsplit, edge_range] : parent_to_child_range_) {
    std::cout << subsplit.SubsplitToString() << ": [" << edge_range.first << ", "
              << edge_range.second << "]" << std::endl;
  }
}

void SubsplitDAG::ToDot(const std::string file_path, bool show_index_labels) const {
  std::ofstream out_file(file_path);
  out_file << ToDot(show_index_labels);
  out_file.close();
}

std::string SubsplitDAG::ToDot(bool show_index_labels) const {
  std::stringstream string_stream;
  string_stream << "digraph g {\n";
  string_stream << "node [shape=record];\n";
  string_stream << "edge [colorscheme=dark23];\n";
  DepthFirstWithAction(
      {GetDAGRootNodeId()},
      SubsplitDAGTraversalAction(
          // BeforeNode
          [this, &string_stream, &show_index_labels](NodeId node_id) {
            const auto &node = GetDAGNode(node_id);
            if (node.IsDAGRootNode()) {
              string_stream << node_id << " [label=\"<f0>&rho;\"]\n";
              return;
            }
            auto bs = node.GetBitset();
            string_stream
                << node_id << " [label=\"<f0>"
                << bs.SubsplitGetClade(SubsplitClade::Left).ToVectorOfSetBitsAsString()
                << "|<f1>";
            if (show_index_labels) {
              string_stream << node_id;
            }
            string_stream
                << "|<f2>"
                << bs.SubsplitGetClade(SubsplitClade::Right).ToVectorOfSetBitsAsString()
                << "\"]\n";
          },
          // AfterNode
          [](NodeId node_id) {},
          // BeforeNodeClade
          [](NodeId node_id, bool is_edge_on_left) {},
          // VisitEdge
          [this, &string_stream, &show_index_labels](NodeId node_id, NodeId child_id,
                                                     bool is_edge_on_left) {
            if (GetDAGNode(child_id).IsLeaf()) {
              string_stream << child_id << " [label=\"<f1>" << child_id << "\"]\n";
            }
            string_stream << "\"" << node_id << "\":";
            string_stream << (is_edge_on_left ? "f0" : "f2");
            string_stream << "->\"";
            string_stream << child_id << "\":f1";
            if (show_index_labels) {
              string_stream << " [label=\"" << GetEdgeIdx(node_id, child_id);
              if (is_edge_on_left) {
                string_stream << "\", color=1, fontcolor=1";
              } else {
                string_stream << "\", color=3, fontcolor=3";
              }
              if (GetDAGNode(node_id).IsDAGRootNode()) {
                string_stream << ",style=dashed]";
              } else {
                string_stream << "]";
              }
            } else {
              if (GetDAGNode(node_id).IsDAGRootNode()) {
                string_stream << "[style=dashed]";
              }
            }
            string_stream << "\n";
          }));
  string_stream << "}";
  return string_stream.str();
}

// ** Build output indexes/vectors methods:

BitsetSizeMap SubsplitDAG::BuildEdgeIndexer() const {
  auto edge_indexer = BitsetSizeMap();
  TopologicalEdgeTraversal([this, &edge_indexer](NodeId parent_id, bool is_edge_on_left,
                                                 NodeId child_id, EdgeId edge_idx) {
    const auto parent_subsplit = GetDAGNode(parent_id).GetBitset(is_edge_on_left);
    const auto child_subsplit = GetDAGNode(child_id).GetBitset();
    SafeInsert(edge_indexer, Bitset::PCSP(parent_subsplit, child_subsplit),
               edge_idx.value_);
  });
  return edge_indexer;
}

SizeBitsetMap SubsplitDAG::BuildInverseEdgeIndexer() const {
  auto edge_to_pcsp_map = SizeBitsetMap();
  TopologicalEdgeTraversal([this, &edge_to_pcsp_map](NodeId parent_id,
                                                     bool is_edge_on_left,
                                                     NodeId child_id, EdgeId edge_idx) {
    const auto parent_subsplit = GetDAGNode(parent_id).GetBitset(is_edge_on_left);
    const auto child_subsplit = GetDAGNode(child_id).GetBitset();
    SafeInsert(edge_to_pcsp_map, edge_idx.value_,
               Bitset::PCSP(parent_subsplit, child_subsplit));
  });
  return edge_to_pcsp_map;
}

TaxonId SubsplitDAG::GetTaxonId(const std::string &name) const {
  Assert(
      ContainsTaxon(name),
      "SubsplitDAG::GetTaxonId(): Taxon with given name is not contained in the DAG.");
  return TaxonId(dag_taxa_.find(name)->second);
}

SubsplitDAGNode SubsplitDAG::GetDAGNode(const NodeId node_id) const {
  Assert(ContainsNode(node_id), "Node with the given node_id does not exist in DAG.");
  return storage_.GetVertices().at(node_id.value_);
}

MutableSubsplitDAGNode SubsplitDAG::GetDAGNode(const NodeId node_id) {
  Assert(ContainsNode(node_id), "Node with the given node_id does not exist in DAG.");
  return storage_.GetVertices().at(node_id.value_);
}

Bitset SubsplitDAG::GetDAGNodeBitset(const NodeId node_id) const {
  return GetDAGNode(node_id).GetBitset();
}

NodeId SubsplitDAG::GetDAGNodeId(const Bitset &subsplit) const {
  Assert(ContainsNode(subsplit), "Node with the given subsplit does not exist in DAG.");
  if (storage_.HaveHost()) {
    auto node = storage_.FindVertex(subsplit);
    return NodeId(node->get().GetId());
  }
  return subsplit_to_id_.at(subsplit);
}

NodeId SubsplitDAG::GetDAGRootNodeId() const { return NodeId(NodeCount() - 1); }

ConstNeighborsView SubsplitDAG::GetRootsplitNodeIds() const {
  return GetDAGNode(GetDAGRootNodeId()).GetLeftLeafward();
}

NodeIdVector SubsplitDAG::GetLeafNodeIds() const {
  NodeIdVector leaf_node_ids;
  for (NodeId node_id = 0; node_id < TaxonCount(); node_id++) {
    leaf_node_ids.push_back(node_id);
  }
  return leaf_node_ids;
}

ConstLineView SubsplitDAG::GetDAGEdge(const EdgeId edge_id) const {
  Assert(ContainsEdge(edge_id), "Node with the given node_id does not exist in DAG.");
  return storage_.GetLine(edge_id).value();
}

Bitset SubsplitDAG::GetDAGEdgeBitset(const EdgeId edge_id) const {
  auto edge = GetDAGEdge(edge_id);
  Bitset parent = GetDAGNodeBitset(edge.GetParent());
  Bitset child = GetDAGNodeBitset(edge.GetChild());
  return Bitset::PCSP(parent, child);
}

EdgeId SubsplitDAG::GetEdgeIdx(const Bitset &parent_subsplit,
                               const Bitset &child_subsplit) const {
  return GetEdgeIdx(GetDAGNodeId(parent_subsplit), GetDAGNodeId(child_subsplit));
}

EdgeId SubsplitDAG::GetEdgeIdx(const NodeId parent_id, const NodeId child_id) const {
  auto edge = storage_.GetLine(parent_id, child_id);
  if (!edge.has_value()) Failwith("Edge not found");
  return EdgeId(edge.value().GetId());
}

EdgeId SubsplitDAG::GetEdgeIdx(const Bitset &edge_pcsp) const {
  return GetEdgeIdx(edge_pcsp.PCSPGetParentSubsplit(),
                    edge_pcsp.PCSPGetChildSubsplit());
}

EdgeIdPair SubsplitDAG::GetChildEdgeRange(const Bitset &subsplit,
                                          const bool is_edge_on_left) const {
  Assert(
      ContainsNode(subsplit),
      "Node with the given subsplit does not exist in SubsplitDAG::GetChildEdgeRange.");
  return parent_to_child_range_.at(SubsplitToSortedOrder(subsplit, is_edge_on_left));
}

StringVector SubsplitDAG::GetSortedVectorOfTaxonNames() const {
  StringVector taxa;
  for (const auto &name_id : dag_taxa_) {
    taxa.push_back(name_id.first);
  }
  std::sort(taxa.begin(), taxa.end());
  return taxa;
}

std::vector<Bitset> SubsplitDAG::GetSortedVectorOfNodeBitsets() const {
  std::vector<Bitset> nodes;
  for (NodeId node_id = NodeId(0); node_id < NodeCount(); node_id++) {
    Bitset node_bitset = GetDAGNode(node_id).GetBitset();
    nodes.push_back(node_bitset);
  }
  std::sort(nodes.begin(), nodes.end());
  return nodes;
}

std::vector<Bitset> SubsplitDAG::GetSortedVectorOfEdgeBitsets() const {
  std::vector<Bitset> edges;
  for (auto i : storage_.GetLines()) {
    auto parent_bitset = GetDAGNode(NodeId(i.GetParent())).GetBitset();
    auto child_bitset = GetDAGNode(NodeId(i.GetChild())).GetBitset();
    Bitset edge_bitset = Bitset::PCSP(parent_bitset, child_bitset);
    edges.push_back(edge_bitset);
  }
  std::sort(edges.begin(), edges.end());
  return edges;
}

const StringTaxonIdMap &SubsplitDAG::GetTaxonMap() const { return dag_taxa_; }

const BitsetNodeIdMap &SubsplitDAG::GetSubsplitToIdMap() const {
  return subsplit_to_id_;
}

EigenVectorXd SubsplitDAG::BuildUniformOnTopologicalSupportPrior() const {
  EigenVectorXd q = EigenVectorXd::Ones(EdgeCountWithLeafSubsplits());

  for (const auto &node_id : RootwardNodeTraversalTrace(true)) {
    const auto &node = GetDAGNode(node_id);
    for (const bool rotated : {false, true}) {
      if (!node.GetLeafward(rotated).empty()) {
        double per_rotated_count = 0.;
        for (auto child_id : node.GetLeafward(rotated)) {
          per_rotated_count += topology_count_below_[child_id.value_];
        }
        for (auto child_id : node.GetLeafward(rotated)) {
          auto edge_idx = GetEdgeIdx(node.Id(), NodeId(child_id));
          q(size_t(edge_idx)) =
              topology_count_below_(child_id.value_) / per_rotated_count;
        }
      }
    }
  }
  return q;
}

Node::NodePtrVec SubsplitDAG::GenerateAllTopologies() const {
  std::vector<Node::NodePtrVec> topology_below(NodeCount());

  auto GetSubtopologies = [&topology_below](SubsplitDAGNode node) {
    Node::NodePtrVec rotated_subtopologies, sorted_subtopologies;
    for (const bool rotated : {false, true}) {
      for (const auto &child_id : node.GetLeafward(rotated)) {
        for (const auto &subtopology : topology_below.at(child_id.value_)) {
          rotated ? rotated_subtopologies.push_back(subtopology)
                  : sorted_subtopologies.push_back(subtopology);
        }
      }
    }
    return std::make_pair(rotated_subtopologies, sorted_subtopologies);
  };

  auto MergeTopologies = [this](NodeId node_id, Node::NodePtrVec &rotated_subtopologies,
                                Node::NodePtrVec &sorted_subtopologies) {
    Node::NodePtrVec topologies;
    for (const auto &rotated_subtopology : rotated_subtopologies) {
      for (const auto &sorted_subtopology : sorted_subtopologies) {
        Node::NodePtr new_topology =
            Node::Join(sorted_subtopology, rotated_subtopology, node_id.value_);
        topologies.push_back(new_topology);
      }
    }
    if (node_id == GetDAGRootNodeId()) {
      // DAG root node has no `sorted_subtopologies`, so loop above yields empty
      // `topologies` vector.
      return rotated_subtopologies;
    }
    return topologies;
  };

  for (const auto &node_id : RootwardNodeTraversalTrace(true)) {
    const auto &node = GetDAGNode(node_id);
    if (node.IsLeaf()) {
      topology_below.at(node_id.value_).push_back(Node::Leaf(node_id.value_));
    } else {
      auto [rotated_topologies, sorted_topologies] = GetSubtopologies(node);
      topology_below[node_id.value_] =
          MergeTopologies(node_id, rotated_topologies, sorted_topologies);
    }
  }

  const auto &topologies = topology_below.at(GetDAGRootNodeId().value_);
  Assert(topologies.size() == TopologyCount(),
         "The realized number of topologies does not match the expected count.");

  // We return a deep copy of every Polished topology to avoid loops in the pointer
  // structure. Such loops can create problems when we Polish the topologies one at a
  // time: polishing a second topology can change the numbering of a} previous topology.
  // This is checked for in the "GPInstance: GenerateCompleteRootedTreeCollection" test.
  Node::NodePtrVec final_topologies;
  final_topologies.reserve(topologies.size());
  for (auto &topology : topologies) {
    topology->Polish();
    final_topologies.push_back(topology->DeepCopy());
  }

  return final_topologies;
}

EigenVectorXd SubsplitDAG::BuildUniformOnAllTopologiesPrior() const {
  EigenVectorXd result = EigenVectorXd::Zero(EdgeCountWithLeafSubsplits());
  for (const auto &[parent_child_id, edge_idx] : storage_.GetLines()) {
    const auto &[parent_id, child_id] = parent_child_id;
    std::ignore = parent_id;
    // If child is a leaf and subsplit is sorted, then child0 will have a zero taxon
    // count.
    auto child_left_taxon_count =
        GetDAGNode(child_id).GetBitset().SubsplitGetClade(SubsplitClade::Left).Count();
    // As long as subsplit is sorted and nonempty, then child1 will have a nonzero taxon
    // count.
    auto child_right_taxon_count =
        GetDAGNode(child_id).GetBitset().SubsplitGetClade(SubsplitClade::Right).Count();
    // The ordering of this subsplit is flipped so that this ratio will be nonzero in
    // the denominator in the case of root & leaves.
    result(size_t(edge_idx)) = Combinatorics::LogChildSubsplitCountRatio(
        child_right_taxon_count, child_left_taxon_count);
  }
  NumericalUtils::Exponentiate(result);

  return result;
}

// ** DAG Lambda Iterator methods:

void SubsplitDAG::IterateOverRealNodes(const NodeLambda &f) const {
  Assert(taxon_count_ < NodeCount(), "No real DAG nodes!");
  for (auto it = storage_.GetVertices().cbegin() + taxon_count_;
       it < storage_.GetVertices().cend() - 1; it++) {
    f(*it);
  }
}

void SubsplitDAG::IterateOverLeafwardEdges(SubsplitDAGNode node, bool rotated,
                                           const NodeLambda &f) const {
  for (const auto child_id : node.GetLeafward(rotated)) {
    f(GetDAGNode(child_id));
  }
}

void SubsplitDAG::IterateOverLeafwardEdges(SubsplitDAGNode node,
                                           const EdgeDestinationLambda &f) const {
  for (bool is_edge_on_left : {false, true}) {
    for (const auto child_id : node.GetLeafward(is_edge_on_left)) {
      f(is_edge_on_left, GetDAGNode(child_id));
    }
  }
}

void SubsplitDAG::IterateOverLeafwardEdgesAndChildren(
    SubsplitDAGNode node, const EdgeAndNodeLambda &f) const {
  IterateOverLeafwardEdges(
      node, [this, &node, &f](bool is_edge_on_left, SubsplitDAGNode child) {
        f(GetEdgeIdx(node.Id(), child.Id()), is_edge_on_left, child.Id());
      });
}

void SubsplitDAG::IterateOverRootwardEdges(SubsplitDAGNode node,
                                           const EdgeDestinationLambda &f) const {
  if (!node.IsRootsplit()) {
    for (bool is_edge_on_left : {false, true}) {
      for (const auto parent_id : node.GetRootward(is_edge_on_left)) {
        f(is_edge_on_left, GetDAGNode(parent_id));
      }
    }
  }
}

void SubsplitDAG::IterateOverRootwardEdgesAndParents(SubsplitDAGNode node,
                                                     const EdgeAndNodeLambda &f) const {
  IterateOverRootwardEdges(
      node, [this, &node, &f](bool rotated, SubsplitDAGNode parent) {
        f(GetEdgeIdx(parent.Id(), node.Id()), rotated, parent.Id());
      });
}

void SubsplitDAG::IterateOverParentAndChildAndLeafwardEdges(
    SubsplitDAGNode node, const ParentRotationChildEdgeLambda &f) const {
  IterateOverLeafwardEdges(
      node, [this, &node, &f](bool is_edge_on_left, SubsplitDAGNode child) {
        f(node.Id(), is_edge_on_left, child.Id(), GetEdgeIdx(node.Id(), child.Id()));
      });
}

RootedIndexerRepresentation SubsplitDAG::IndexerRepresentationOf(
    const BitsetSizeMap &indexer, const Node::NodePtr &topology,
    size_t default_index) const {
  return RootedSBNMaps::IndexerRepresentationOf(indexer, topology, default_index);
}

EigenVectorXd SubsplitDAG::UnconditionalNodeProbabilities(
    EigenConstVectorXdRef normalized_sbn_parameters) const {
  EigenVectorXd node_probabilities(NodeCount());
  node_probabilities.setZero();
  node_probabilities[GetDAGRootNodeId().value_] = 1.;

  TopologicalEdgeTraversal([&node_probabilities, &normalized_sbn_parameters](
                               const NodeId parent_id, const bool is_edge_on_left,
                               const NodeId child_id, const EdgeId edge_idx) {
    const double child_probability_given_parent =
        normalized_sbn_parameters[edge_idx.value_];
    Assert(child_probability_given_parent >= 0. && child_probability_given_parent <= 1.,
           "UnconditionalNodeProbabilities: got an out-of-range probability. Are these "
           "normalized and in linear space?");
    const double parent_probability = node_probabilities[parent_id.value_];
    node_probabilities[child_id.value_] +=
        parent_probability * child_probability_given_parent;
  });

  return node_probabilities;
}

BitsetDoubleMap SubsplitDAG::UnconditionalSubsplitProbabilities(
    EigenConstVectorXdRef normalized_sbn_parameters) const {
  auto node_probabilities = UnconditionalNodeProbabilities(normalized_sbn_parameters);
  BitsetDoubleMap subsplit_probability_map;
  for (NodeId node_id = NodeId(0);
       static_cast<Eigen::Index>(node_id.value_) < node_probabilities.size();
       node_id++) {
    const auto &subsplit_bitset = GetDAGNode(node_id).GetBitset();
    if (node_id != GetDAGRootNodeId() && !subsplit_bitset.SubsplitIsLeaf()) {
      SafeInsert(subsplit_probability_map, subsplit_bitset,
                 node_probabilities[node_id.value_]);
    }
  }
  return subsplit_probability_map;
}

EigenVectorXd SubsplitDAG::InvertedGPCSPProbabilities(
    EigenConstVectorXdRef normalized_sbn_parameters,
    EigenConstVectorXdRef node_probabilities) const {
  EigenVectorXd inverted_probabilities =
      EigenVectorXd(normalized_sbn_parameters.size());
  inverted_probabilities.setOnes();
  TopologicalEdgeTraversal(
      [this, &node_probabilities, &normalized_sbn_parameters, &inverted_probabilities](
          const NodeId parent_id, const bool is_edge_on_left, const NodeId child_id,
          const EdgeId edge_idx) {
        // The traversal doesn't set the rootsplit probabilities, but those are always 1
        // (there is only one "parent" of a rootsplit).
        if (parent_id != GetDAGRootNodeId()) {
          // For a PCSP t -> s:
          inverted_probabilities[edge_idx.value_] =         // P(t|s)
              node_probabilities[parent_id.value_] *        // P(t)
              normalized_sbn_parameters[edge_idx.value_] /  // P(s|t)
              node_probabilities[child_id.value_];          // P(s)
        }
      });
  return inverted_probabilities;
}

BitsetVector SubsplitDAG::GetChildSubsplits(const SizeBitsetMap &index_to_child,
                                            const Bitset &parent_subsplit,
                                            bool include_leaf_subsplits) {
  BitsetVector children_subsplits;
  // Add all non-leaf child subsplit bitsets.
  if (parent_to_child_range_.count(parent_subsplit) > 0) {
    const auto [start, stop] = parent_to_child_range_.at(parent_subsplit);
    for (auto idx = start; idx < stop; idx++) {
      children_subsplits.push_back(index_to_child.at(idx.value_));
    }
  }
  // Optionally add leaf child subsplit bitsets.
  else if (include_leaf_subsplits) {
    // This method is designed to be called before calling
    // AddLeafSubsplitsToDAGEdgesAndParentToRange. In that case, if the second clade
    // of the subsplit is just a single taxon, the subsplit will not map to any value in
    // parent_to_child_range_.
    //
    // But we still need to create and connect to leaf subsplits in the DAG. So, here we
    // make a leaf child subsplit.
    children_subsplits.push_back(Bitset::LeafSubsplitOfParentSubsplit(parent_subsplit));
  }

  return children_subsplits;
}

// ** Modify DAG methods:

BitsetEdgeIdPairMap BitsetSizePairMapToBitsetEdgeIdPairMap(
    BitsetSizePairMap &bitset_size_map) {
  BitsetEdgeIdPairMap bitset_edgeid_map;
  for (const auto &[bitset, id_pair] : bitset_size_map) {
    const auto &[begin_id, end_id] = id_pair;
    bitset_edgeid_map[bitset] = {EdgeId(begin_id), EdgeId(end_id)};
  }
  return bitset_edgeid_map;
}

std::tuple<BitsetSizeMap, SizeBitsetMap, BitsetVector>
SubsplitDAG::ProcessTopologyCounter(const Node::TopologyCounter &topology_counter) {
  BitsetSizeMap edge_indexer;
  SizeBitsetMap index_to_child;
  BitsetVector rootsplits;
  BitsetSizePairMap parent_to_child_range;
  std::tie(rootsplits, edge_indexer, index_to_child, parent_to_child_range,
           edge_count_without_leaf_subsplits_) =
      SBNMaps::BuildIndexerBundle(RootedSBNMaps::RootsplitCounterOf(topology_counter),
                                  RootedSBNMaps::PCSPCounterOf(topology_counter));
  parent_to_child_range_ =
      BitsetSizePairMapToBitsetEdgeIdPairMap(parent_to_child_range);
  return {edge_indexer, index_to_child, rootsplits};
}

void SubsplitDAG::BuildTaxonMap(const TagStringMap &tag_taxon_map) {
  // Insert all taxa from tree_collections's map to SubsplitDAG map.
  for (const auto &[tag, name] : tag_taxon_map) {
    // The "tag" key of the tree_collection's taxon_map is 2 bitpacked ints: [id,
    // topology count]. We only care about the id.
    TaxonId id = TaxonId{static_cast<size_t>(UnpackFirstInt(tag))};
    dag_taxa_.insert(std::make_pair(name, id));
  }
}

NodeId SubsplitDAG::CreateAndInsertNode(const Bitset &subsplit) {
  NodeId node_id = NodeId(NodeCount());
  storage_.AddVertex({node_id, subsplit});
  SafeInsert(subsplit_to_id_, subsplit, node_id);
  return node_id;
}

EdgeId SubsplitDAG::CreateAndInsertEdge(const NodeId parent_id, const NodeId child_id,
                                        const bool is_edge_on_left) {
  Assert(ContainsNode(parent_id), "Node with the given parent_id does not exist.");
  Assert(ContainsNode(child_id), "Node with the given child_id does not exist.");
  // Insert edge between parent and child.
  EdgeId edge_idx = EdgeId(EdgeCountWithLeafSubsplits());
  ConnectGivenNodes(parent_id, child_id, is_edge_on_left, edge_idx);
  storage_.AddLine({edge_idx, parent_id, child_id,
                    is_edge_on_left ? SubsplitClade::Left : SubsplitClade::Right});
  return edge_idx;
}

void SubsplitDAG::ConnectGivenNodes(const NodeId parent_id, const NodeId child_id,
                                    const bool is_edge_on_left, const EdgeId edge_id) {
  auto parent_node = GetDAGNode(parent_id);
  auto child_node = GetDAGNode(child_id);

  SubsplitClade which_parent_clade =
      (is_edge_on_left ? SubsplitClade::Left : SubsplitClade::Right);
  parent_node.AddEdge(child_node.Id(), edge_id, Direction::Leafward,
                      which_parent_clade);
  child_node.AddEdge(parent_node.Id(), edge_id, Direction::Rootward,
                     which_parent_clade);
}

void SubsplitDAG::ConnectNodes(const SizeBitsetMap &index_to_child, NodeId node_id,
                               bool is_edge_on_left) {
  // Get bitset of parent node according to its rotation.
  const auto subsplit = GetDAGNode(node_id).GetBitset(is_edge_on_left);
  // Build vector of child node's subsplits.
  const auto children = GetChildSubsplits(index_to_child, subsplit, true);
  // Connect parent node to all child nodes.
  for (const auto &child_subsplit : children) {
    ConnectGivenNodes(node_id, GetDAGNodeId(child_subsplit), is_edge_on_left,
                      EdgeId(NoId));
  }
}

void SubsplitDAG::BuildNodes(const SizeBitsetMap &index_to_child,
                             const BitsetVector &rootsplits) {
  std::unordered_set<Bitset> visited_subsplits;
  // We will create leaf subsplits and insert to storage_ nodes.
  // Leaf subsplits (i.e. leaf nodes) will take IDs in [0, taxon_count_).
  for (size_t taxon_idx = 0; taxon_idx < taxon_count_; taxon_idx++) {
    CreateAndInsertNode(Bitset::LeafSubsplitOfNonemptyClade(
        Bitset::Singleton(taxon_count_, taxon_idx)));
  }
  // Then we add the remaining nodes.
  // The root splits will take on the higher IDs compared to the non-rootsplits.
  for (const auto &rootsplit : rootsplits) {
    BuildNodesDepthFirst(index_to_child, rootsplit, visited_subsplits);
  }
  // #350 look for "DAG root" and fix
  // Finally, we add the DAG root node.
  CreateAndInsertNode(Bitset::UCASubsplitOfTaxonCount(taxon_count_));
}

void SubsplitDAG::BuildNodesDepthFirst(const SizeBitsetMap &index_to_child,
                                       const Bitset &subsplit,
                                       std::unordered_set<Bitset> &visited_subsplits) {
  visited_subsplits.insert(subsplit);
  for (bool rotated : {false, true}) {
    for (const auto &child_subsplit : GetChildSubsplits(
             index_to_child, SubsplitToSortedOrder(subsplit, rotated), false)) {
      if (visited_subsplits.count(child_subsplit) == 0) {
        BuildNodesDepthFirst(index_to_child, child_subsplit, visited_subsplits);
      }
    }
  }
  CreateAndInsertNode(subsplit);
}

void SubsplitDAG::BuildEdges(const SizeBitsetMap &index_to_child) {
  // Connect every node except for the DAG root node.
  for (NodeId node_id = NodeId(taxon_count_); node_id < GetDAGRootNodeId(); node_id++) {
    ConnectNodes(index_to_child, node_id, true);
    ConnectNodes(index_to_child, node_id, false);
  }
  // Connect the DAG root node.
  ConnectNodes(index_to_child, GetDAGRootNodeId(), true);
}

void SubsplitDAG::BuildDAGEdgesFromEdgeIndexer(BitsetSizeMap &edge_indexer) {
  for (const auto &[edge, index] : edge_indexer) {
    Assert(edge.size() == 3 * taxon_count_,
           "All edges should be bitsets with size 3 times taxon_count_.");
    const auto parent_id = GetDAGNodeId(edge.PCSPGetParentSubsplit());
    const auto child_id = GetDAGNodeId(edge.PCSPGetChildSubsplit());
    auto left = GetDAGNode(parent_id).GetLeftLeafward();
    auto clade = std::find(left.begin(), left.end(), child_id) != left.end()
                     ? SubsplitClade::Left
                     : SubsplitClade::Right;
    if (clade == SubsplitClade::Right) {
      auto right = GetDAGNode(parent_id).GetRightLeafward();
      Assert(std::find(right.begin(), right.end(), child_id) != right.end(),
             "Parent has no connection to child");
    }
    storage_.AddLine({EdgeId(index), parent_id, child_id, clade});
  }
}

void SubsplitDAG::AddLeafSubsplitsToDAGEdgesAndParentToRange() {
  for (NodeId node_id = NodeId(0); node_id < taxon_count_; node_id++) {
    const auto current_bitset(GetDAGNode(node_id).GetBitset());
    IterateOverRootwardEdges(GetDAGNode(node_id), [this, current_bitset](
                                                      const bool is_edge_on_left,
                                                      SubsplitDAGNode node) {
      SafeInsert(parent_to_child_range_, node.GetBitset(is_edge_on_left),
                 {EdgeId(EdgeCountWithLeafSubsplits()),
                  EdgeId(EdgeCountWithLeafSubsplits() + 1)});
      storage_.AddLine({EdgeId(EdgeCountWithLeafSubsplits()), node.Id(),
                        GetDAGNodeId(current_bitset),
                        is_edge_on_left ? SubsplitClade::Left : SubsplitClade::Right});
    });
  }
}

void SubsplitDAG::StoreEdgeIds() {
  for (auto edge : storage_.GetLines()) {
    auto parent = storage_.GetVertices().at(edge.GetParent().value_);
    auto child = storage_.GetVertices().at(edge.GetChild().value_);
    parent.SetEdgeId(edge.GetChild(), edge.GetId());
    child.SetEdgeId(edge.GetParent(), edge.GetId());
  }
}

void SubsplitDAG::RootwardDepthFirst(NodeId node_id, NodeIdVector &visit_order,
                                     std::unordered_set<NodeId> &visited_nodes) const {
  // Add to set of all visited nodes.
  SafeInsert(visited_nodes, node_id);
  // #350: Look for sorted/rotated and update.
  // Recurse on sorted children.
  const auto &node = GetDAGNode(node_id);
  for (auto parent_id : node.GetRightRootward()) {
    if (visited_nodes.count(parent_id) == 0) {
      RootwardDepthFirst(parent_id, visit_order, visited_nodes);
    }
  }
  // Recurse on rotated children.
  for (auto parent_id : node.GetLeftRootward()) {
    if (visited_nodes.count(parent_id) == 0) {
      RootwardDepthFirst(parent_id, visit_order, visited_nodes);
    }
  }
  // Append to vector post-order (after all children have been visited).
  visit_order.push_back(node_id);
}

void SubsplitDAG::LeafwardDepthFirst(NodeId node_id, NodeIdVector &visit_order,
                                     std::unordered_set<NodeId> &visited_nodes) const {
  // Add to set of all visited nodes.
  SafeInsert(visited_nodes, node_id);
  // Recurse on right/sorted children.
  for (auto child_id : GetDAGNode(node_id).GetRightLeafward()) {
    if (visited_nodes.count(child_id) == 0) {
      LeafwardDepthFirst(child_id, visit_order, visited_nodes);
    }
  }
  // Recurse on left/rotated children.
  for (auto child_id : GetDAGNode(node_id).GetLeftLeafward()) {
    if (visited_nodes.count(child_id) == 0) {
      LeafwardDepthFirst(child_id, visit_order, visited_nodes);
    }
  }
  // Append to vector post-order (after all children have been visited).
  visit_order.push_back(node_id);
}

NodeIdVector SubsplitDAG::LeafwardNodeTraversalTrace(bool include_dag_root_node) const {
  NodeIdVector visit_order;
  std::unordered_set<NodeId> visited_nodes;
  if (!include_dag_root_node) {
    SafeInsert(visited_nodes, GetDAGRootNodeId());
  }
  for (NodeId leaf_id = NodeId(0); leaf_id < taxon_count_; leaf_id++) {
    RootwardDepthFirst(leaf_id, visit_order, visited_nodes);
  }
  return visit_order;
}

NodeIdVector SubsplitDAG::RootwardNodeTraversalTrace(bool include_dag_root_node) const {
  NodeIdVector visit_order;
  std::unordered_set<NodeId> visited_nodes;
  for (const auto &rootsplit_id : GetRootsplitNodeIds()) {
    LeafwardDepthFirst(rootsplit_id, visit_order, visited_nodes);
  }
  if (include_dag_root_node) {
    visit_order.push_back(GetDAGRootNodeId());
  }
  return visit_order;
}

NodeIdVector SubsplitDAG::TopologicalNodeTraversalTrace() const {
  auto visit_order = RootwardNodeTraversalTrace(true);
  std::reverse(visit_order.begin(), visit_order.end());
  return visit_order;
}

void SubsplitDAG::TopologicalEdgeTraversal(ParentRotationChildEdgeLambda f) const {
  for (const auto node_id : TopologicalNodeTraversalTrace()) {
    IterateOverLeafwardEdgesAndChildren(
        GetDAGNode(node_id),
        [&f, &node_id](const EdgeId edge_idx, const bool is_edge_on_left,
                       const NodeId child_id) {
          f(node_id, is_edge_on_left, child_id, edge_idx);
        });
  }
}

// ** Miscellaneous

Bitset SubsplitDAG::SubsplitToSortedOrder(const Bitset &subsplit, bool rotated) const {
  return rotated ? subsplit.SubsplitRotate() : subsplit;
}

SizeVector SubsplitDAG::BuildTaxonTranslationMap(const SubsplitDAG &dag_a,
                                                 const SubsplitDAG &dag_b) {
  auto names_a = dag_a.GetSortedVectorOfTaxonNames();
  auto names_b = dag_b.GetSortedVectorOfTaxonNames();
  Assert(names_a == names_b,
         "SubsplitDAG::BuildTaxonTranslationMap(): SubsplitDAGs do not cover the same "
         "taxon set.");
  Assert(names_a.size() == dag_a.TaxonCount(),
         "SubsplitDAG::BuildTaxonTranslationMap(): Number of taxon names does not "
         "match the number of taxa in the DAG.");

  SizeVector taxon_map(names_a.size());
  for (const auto &name : names_a) {
    taxon_map[dag_a.GetTaxonId(name).value_] = dag_b.GetTaxonId(name).value_;
  }
  return taxon_map;
}

int SubsplitDAG::BitsetCompareViaTaxonTranslationMap(const Bitset &bitset_a,
                                                     const Bitset &bitset_b,
                                                     const SizeVector &taxon_map) {
  Bitset bitset_a_translated_to_b =
      BitsetTranslateViaTaxonTranslationMap(bitset_a, taxon_map);
  return Bitset::Compare(bitset_a_translated_to_b, bitset_b);
}

Bitset SubsplitDAG::BitsetTranslateViaTaxonTranslationMap(
    const Bitset &bitset, const SizeVector &taxon_map, const bool forward_translate) {
  Bitset translated_bitset(bitset.size());
  size_t clade_size = taxon_map.size();
  size_t clade_count = bitset.size() / taxon_map.size();
  // Remap each clade individually.
  for (size_t i = 0; i < clade_count; i++) {
    // Remap each bit in clade according to a new position according to the taxon map.
    for (size_t j = 0; j < clade_size; j++) {
      size_t input_offset = (i * clade_size) + j;
      size_t output_offset = (i * clade_size) + taxon_map[j];
      if (forward_translate) {
        translated_bitset.set(input_offset, bitset[output_offset]);
      } else {
        translated_bitset.set(output_offset, bitset[input_offset]);
      }
    }
  }
  return translated_bitset;
}

TagStringMap SubsplitDAG::BuildDummyTagTaxonMap(const size_t taxon_count) {
  TagStringMap tag_taxon_map;
  for (size_t i = 0; i < taxon_count; i++) {
    std::string name = std::string("x" + std::to_string(i));
    tag_taxon_map.insert(std::make_pair(PackInts(i, 0), name));
  }
  return tag_taxon_map;
}

// ** Query DAG

bool SubsplitDAG::ContainsTaxon(const std::string &name) const {
  return dag_taxa_.find(name) != dag_taxa_.end();
}

bool SubsplitDAG::ContainsNode(const Bitset &subsplit) const {
  if (storage_.HaveHost()) {
    if (storage_.FindVertex(subsplit).has_value()) return true;
  }
  return subsplit_to_id_.find(subsplit) != subsplit_to_id_.end();
}

bool SubsplitDAG::ContainsNode(const NodeId node_id) const {
  return node_id < NodeCount();
}

bool SubsplitDAG::ContainsEdge(const Bitset &parent_subsplit,
                               const Bitset &child_subsplit) const {
  if (!(ContainsNode(parent_subsplit) && ContainsNode(child_subsplit))) {
    return false;
  }
  return ContainsEdge(GetDAGNodeId(parent_subsplit), GetDAGNodeId(child_subsplit));
}

bool SubsplitDAG::ContainsEdge(const NodeId parent_id, const NodeId child_id) const {
  return storage_.GetLine(parent_id, child_id).has_value();
}

bool SubsplitDAG::ContainsEdge(const Bitset &edge_pcsp) const {
  return ContainsEdge(edge_pcsp.PCSPGetParentSubsplit(),
                      edge_pcsp.PCSPGetChildSubsplit());
}

bool SubsplitDAG::ContainsEdge(const EdgeId edge_id) const {
  return storage_.GetLine(edge_id).has_value();
}

bool SubsplitDAG::IsNodeRoot(const NodeId node_id) const {
  return (node_id == GetDAGRootNodeId());
}

bool SubsplitDAG::IsNodeLeaf(const NodeId node_id) const {
  return (node_id < TaxonCount());
}

bool SubsplitDAG::IsEdgeRoot(const EdgeId edge_id) const {
  const auto parent_id = NodeId(GetDAGEdge(edge_id).GetParent());
  return IsNodeRoot(parent_id);
}

bool SubsplitDAG::IsEdgeLeaf(const EdgeId edge_id) const {
  const auto child_id = NodeId(GetDAGEdge(edge_id).GetChild());
  return IsNodeLeaf(child_id);
}

bool SubsplitDAG::ContainsTopology(const Node::NodePtr topology,
                                   const bool is_quiet) const {
  std::stringstream dev_null;
  std::ostream &os = (is_quiet ? dev_null : std::cerr);
  bool contains_topology = true;
  BoolVector leaf_check(TaxonCount(), false);
  // iterate over rest of topology.
  topology->Preorder([this, &os, &leaf_check, &contains_topology](const Node *node) {
    // Skip if already found topology not in DAG.
    if (!contains_topology) {
      return;
    }
    // Check that topology covers entire taxon set.
    if (node->Leaves().size() != TaxonCount()) {
      os << "DoesNotContainTopology: Number of topology leaves different size than "
            "taxon set."
         << std::endl;
      contains_topology = false;
      return;
    }
    // If node is a child, make sure it is a singleton and check the leaf bit.
    if (node->IsLeaf()) {
      const auto singleton = node->Leaves().SingletonOption();
      if (!singleton.has_value()) {
        os << "DoesNotContainTopology: Leaf node is not a singleton. -- "
           << node->Leaves() << std::endl;
        contains_topology = false;
        return;
      }
      leaf_check[singleton.value()] = true;
    }
    // Otherwise, find both child PCSPs from node and check that they are in the DAG.
    else {
      const auto child_nodes = node->Children();
      if (child_nodes.size() != 2) {
        os << "DoesNotContainTopology: Non-leaf node does not have 2 children."
           << std::endl;
        contains_topology = false;
        return;
      }
      const auto parent_subsplit =
          Bitset::Subsplit(child_nodes[0]->Leaves(), child_nodes[1]->Leaves());
      for (const auto &child_node : child_nodes) {
        const auto grandchild_nodes = child_node->Children();
        const auto child_subsplit =
            (child_node->IsLeaf())
                ? Bitset::Subsplit(child_node->Leaves(), Bitset(TaxonCount(), false))
                : Bitset::Subsplit(grandchild_nodes[0]->Leaves(),
                                   grandchild_nodes[1]->Leaves());
        if (!ContainsEdge(parent_subsplit, child_subsplit)) {
          os << "DoesNotContainTopology: Edge in topology not found in DAG -- "
             << Bitset::PCSP(parent_subsplit, child_subsplit).PCSPToString()
             << std::endl;
          contains_topology = false;
          return;
        }
      }
    }
  });
  if (!contains_topology) {
    return false;
  }
  // Check that every leaf node has been visited.
  bool all_leaves = std::all_of(leaf_check.begin(), leaf_check.end(),
                                [](bool all_true) { return all_true; });
  if (!all_leaves) {
    os << "DoesNotContainTopology: Topology does not span every leaf -- " << leaf_check
       << std::endl;
    return false;
  }
  return true;
}

// ** Build Output Indexers/Vectors

NodeIdVectorPair SubsplitDAG::BuildParentIdVectors(const Bitset &subsplit) const {
  NodeIdVector left_parents, right_parents;
  // Linear search for all parents.
  for (const auto &[potential_parent_subsplit, node_id] : subsplit_to_id_) {
    if (subsplit.SubsplitIsLeftChildOf(potential_parent_subsplit)) {
      left_parents.push_back(node_id);
    } else if (subsplit.SubsplitIsRightChildOf(potential_parent_subsplit)) {
      right_parents.push_back(node_id);
    }
  }
  return {left_parents, right_parents};
}

NodeIdVectorPair SubsplitDAG::BuildChildIdVectors(const Bitset &subsplit) const {
  NodeIdVector left_children, right_children;
  // Linear search for all children.
  for (const auto &[potential_child_subsplit, node_id] : subsplit_to_id_) {
    if (potential_child_subsplit.SubsplitIsLeftChildOf(subsplit)) {
      left_children.push_back(node_id);
    } else if (potential_child_subsplit.SubsplitIsRightChildOf(subsplit)) {
      right_children.push_back(node_id);
    }
  }
  return {left_children, right_children};
}

void SubsplitDAG::ConnectChildToAllChildren(const Bitset &child_subsplit,
                                            EdgeIdVector &added_edge_idxs) {
  const auto [left_leafward_of_child, right_leafward_of_child] =
      BuildChildIdVectors(child_subsplit);

  for (const auto &[children_of_child, rotated] :
       std::vector<std::pair<NodeIdVector, bool>>{{left_leafward_of_child, true},
                                                  {right_leafward_of_child, false}}) {
    SafeInsert(parent_to_child_range_, SubsplitToSortedOrder(child_subsplit, rotated),
               {EdgeId(EdgeCountWithLeafSubsplits()),
                EdgeId(EdgeCountWithLeafSubsplits() + children_of_child.size())});

    for (const auto child_of_child_id : children_of_child) {
      const auto new_edge_idx =
          CreateAndInsertEdge(GetDAGNodeId(child_subsplit), child_of_child_id, rotated);
      added_edge_idxs.push_back(new_edge_idx);
    }
  }
}

void SubsplitDAG::ConnectParentToAllChildrenExcept(const Bitset &parent_subsplit,
                                                   const Bitset &child_subsplit,
                                                   EdgeIdVector &added_edge_idxs) {
  const auto [left_leafward_of_parent, right_leafward_of_parent] =
      BuildChildIdVectors(parent_subsplit);

  for (const auto &[children_of_parent, is_edge_on_left] :
       std::vector<std::pair<NodeIdVector, bool>>{{left_leafward_of_parent, true},
                                                  {right_leafward_of_parent, false}}) {
    SafeInsert(parent_to_child_range_,
               SubsplitToSortedOrder(parent_subsplit, is_edge_on_left),
               {EdgeId(EdgeCountWithLeafSubsplits()),
                EdgeId(EdgeCountWithLeafSubsplits() + children_of_parent.size())});

    for (const auto child_of_parent_id : children_of_parent) {
      if (child_of_parent_id != GetDAGNodeId(child_subsplit)) {
        const auto new_edge_idx = CreateAndInsertEdge(
            GetDAGNodeId(parent_subsplit), child_of_parent_id, is_edge_on_left);
        added_edge_idxs.push_back(new_edge_idx);
      }
    }
  }
}

void SubsplitDAG::ConnectChildToAllParentsExcept(const Bitset &parent_subsplit,
                                                 const Bitset &child_subsplit,
                                                 EdgeIdVector &added_edge_idxs) {
  const auto [left_rootward_of_child, right_rootward_of_child] =
      BuildParentIdVectors(child_subsplit);

  for (const auto &[parents_of_child, rotated] :
       std::vector<std::pair<NodeIdVector, bool>>{{left_rootward_of_child, true},
                                                  {right_rootward_of_child, false}}) {
    for (const auto parent_of_child_id : parents_of_child) {
      if (parent_of_child_id != GetDAGNodeId(parent_subsplit)) {
        const auto new_edge_idx = CreateAndInsertEdge(
            parent_of_child_id, GetDAGNodeId(child_subsplit), rotated);
        added_edge_idxs.push_back(new_edge_idx);
      }
    }
  }
}

void SubsplitDAG::ConnectParentToAllParents(const Bitset &parent_subsplit,
                                            EdgeIdVector &added_edge_idxs) {
  const auto [left_rootward_of_parent, right_rootward_of_parent] =
      BuildParentIdVectors(parent_subsplit);

  for (const auto &[parents_of_parent, rotated] :
       std::vector<std::pair<NodeIdVector, bool>>{{left_rootward_of_parent, true},
                                                  {right_rootward_of_parent, false}}) {
    for (const auto parent_of_parent_id : parents_of_parent) {
      const auto new_edge_idx = CreateAndInsertEdge(
          parent_of_parent_id, GetDAGNodeId(parent_subsplit), rotated);
      added_edge_idxs.push_back(new_edge_idx);
    }
  }
}

SubsplitDAG::ModificationResult SubsplitDAG::AddNodePair(const NNIOperation &nni) {
  return AddNodePair(nni.parent_, nni.child_);
}

SubsplitDAG::ModificationResult SubsplitDAG::AddNodePair(const Bitset &parent_subsplit,
                                                         const Bitset &child_subsplit) {
  // Check that node pair will create a valid SubsplitDAG.
  Assert(IsValidAddNodePair(parent_subsplit, child_subsplit),
         "The given pair of nodes is incompatible with DAG in "
         "SubsplitDAG::AddNodePair.");
  // Initialize output vectors.
  NodeIdVector added_node_ids;
  EdgeIdVector added_edge_idxs;
  Reindexer node_reindexer, edge_reindexer;
  // Check if either parent or child don't already exist in the DAG.
  const bool parent_is_new = !ContainsNode(parent_subsplit);
  const bool child_is_new = !ContainsNode(child_subsplit);
  // Soft assert: This allows for parent-child pair to exist in the DAG, but no work
  // is done. If both the parent and child already exist in DAG, return added_node_ids
  // and added_edge_idxs as empty, and node_reindexer and edge_reindexer as identity
  // reindexers.
  if (!parent_is_new && !child_is_new) {
    // Return default reindexers if both nodes already exist.
    node_reindexer = Reindexer::IdentityReindexer(NodeCount());
    edge_reindexer = Reindexer::IdentityReindexer(EdgeCountWithLeafSubsplits());
    return {added_node_ids, added_edge_idxs, node_reindexer, edge_reindexer};
  }

  // Note: `prev_node_count` acts as a place marker. We know what the DAG root node id
  // is
  // (`prev_node_count - 1`).
  const size_t prev_node_count = NodeCount();

  // Add parent/child nodes and connect them to their children
  // If child node is new, add node and connect it to all its children.
  if (child_is_new) {
    CreateAndInsertNode(child_subsplit);
    added_node_ids.push_back(GetDAGNodeId(child_subsplit));
    // Don't reindex these edges.
    ConnectChildToAllChildren(child_subsplit, added_edge_idxs);
  }
  // If parent node is new, add node it to all its children (except )
  if (parent_is_new) {
    CreateAndInsertNode(parent_subsplit);
    added_node_ids.push_back(GetDAGNodeId(parent_subsplit));
    // Don't reindex these edges.
    ConnectParentToAllChildrenExcept(parent_subsplit, child_subsplit, added_edge_idxs);
  }
  // Note: `prev_edge_count` is a marker conveying where we need to start
  // reindexing edge idxs.
  // Edges are only reindexed if the parent node already existed in the DAG
  // (so as to ensure that edges descending from the same node clade have
  // contiguous idxs).
  size_t prev_edge_count = EdgeCountWithLeafSubsplits();
  // Connect the given parent node to the given child node.
  added_edge_idxs.push_back(EdgeId(EdgeCountWithLeafSubsplits()));
  CreateAndInsertEdge(GetDAGNodeId(parent_subsplit), GetDAGNodeId(child_subsplit),
                      child_subsplit.SubsplitIsLeftChildOf(parent_subsplit));
  // Don't reindex the edge between the given parent and child if the parent is new.
  if (parent_is_new) {
    prev_edge_count = EdgeCountWithLeafSubsplits();
  }
  if (child_is_new) {
    // Reindex these edges.
    ConnectChildToAllParentsExcept(parent_subsplit, child_subsplit, added_edge_idxs);
  }
  if (parent_is_new) {
    // Reindex these edges.
    ConnectParentToAllParents(parent_subsplit, added_edge_idxs);
  }
  if (!storage_.HaveHost()) {
    // Create reindexers.
    node_reindexer = BuildNodeReindexer(prev_node_count);
    edge_reindexer = BuildEdgeReindexer(prev_edge_count);
    // Update the ids in added_node_ids and added_edge_idxs according to the
    // reindexers.
    Reindexer::RemapIdVector<NodeId>(added_node_ids, node_reindexer);
    Reindexer::RemapIdVector<EdgeId>(added_edge_idxs, edge_reindexer);
    // Update fields in the Subsplit DAG according to the reindexers.
    RemapNodeIds(node_reindexer);
    RemapEdgeIdxs(edge_reindexer);
    // Recount topologies.
    CountTopologies();
  }

  return {added_node_ids, added_edge_idxs, node_reindexer, edge_reindexer};
}

SubsplitDAG::ModificationResult SubsplitDAG::FullyConnect() {
  NodeIdVector added_node_ids;
  EdgeIdVector added_edge_idxs;
  Reindexer node_reindexer, edge_reindexer;
  node_reindexer = Reindexer::IdentityReindexer(NodeCount());
  size_t prv_edge_count = EdgeCountWithLeafSubsplits();
  // Find potential edges
  for (const auto &node : storage_.GetVertices()) {
    if (node.IsLeaf()) {
      continue;
    }
    const auto [left_children, right_children] = BuildChildIdVectors(node.GetBitset());
    for (const auto &children : {left_children, right_children}) {
      const bool is_on_left = (children == left_children);
      for (const auto child_id : children) {
        if (!ContainsEdge(node.Id(), child_id)) {
          const auto edge_idx = CreateAndInsertEdge(node.Id(), child_id, is_on_left);
          added_edge_idxs.push_back(edge_idx);
        }
      }
    }
  }
  // Create reindexer and update fields.
  edge_reindexer = BuildEdgeReindexer(prv_edge_count);
  Reindexer::RemapIdVector<EdgeId>(added_edge_idxs, edge_reindexer);
  RemapEdgeIdxs(edge_reindexer);
  // Recount topologies.
  CountTopologies();

  return {added_node_ids, added_edge_idxs, node_reindexer, edge_reindexer};
}

// ** Validation Tests

bool SubsplitDAG::IsConsistent() const {
  Failwith("SubsplitDAG::IsConsistent() is not yet implemented.");
  return false;
}

bool SubsplitDAG::IsValid() const {
  size_t correct_id = 0;
  for (auto node : storage_.GetVertices()) {
    if (correct_id++ != node.Id().value_) {
      return false;
    }
    if (!node.IsValid()) {
      return false;
    }
  }
  return true;
}

bool SubsplitDAG::IsValidAddNodePair(const Bitset &parent_subsplit,
                                     const Bitset &child_subsplit) const {
  // Get the number of adjacent nodes in the given direction.
  auto GetNodeCounts = [this](const Bitset &subsplit, const bool is_rootward) {
    const auto [left, right] =
        is_rootward ? BuildParentIdVectors(subsplit) : BuildChildIdVectors(subsplit);
    SizePair sides = {left.size(), right.size()};
    return sides;
  };
  // Get all adjacent nodes, not including parent and child.
  auto [left_rootward_of_parent, right_rootward_of_parent] =
      GetNodeCounts(parent_subsplit, true);
  auto [left_leafward_of_parent, right_leafward_of_parent] =
      GetNodeCounts(parent_subsplit, false);
  auto [left_leafward_of_child, right_leafward_of_child] =
      GetNodeCounts(child_subsplit, false);
  // Add child to parent's adjacent nodes.
  const bool is_left_child = child_subsplit.SubsplitIsLeftChildOf(parent_subsplit);
  if (is_left_child) {
    left_leafward_of_parent++;
  } else {
    right_leafward_of_parent++;
  }

  // (1) Added nodes are parent/child pair.
  if (Bitset::SubsplitIsParentChildPair(parent_subsplit, child_subsplit) == false) {
    return false;
  }
  // (2) Nodes do not add or remove taxa.
  if ((parent_subsplit.size() != 2 * taxon_count_) ||
      (child_subsplit.size() != 2 * taxon_count_)) {
    return false;
  }
  // (3) The parent node has at least one parent, and at least one left and right
  // child (including the added child node).
  const bool parent_has_parent =
      (left_rootward_of_parent > 0) || (right_rootward_of_parent > 0);
  const bool parent_has_children =
      (left_leafward_of_parent > 0) && (right_leafward_of_parent > 0);
  if ((parent_has_parent && parent_has_children) == false) {
    return false;
  }
  // (4) The child node has at least one parent, and at least one rotated and sorted
  // child.
  // (* we know child node has a parent node, so only need to check children)
  const bool child_has_children =
      (left_leafward_of_child > 0) && (right_leafward_of_child > 0);
  if (child_has_children == false) {
    return false;
  }
  return true;
}

bool SubsplitDAG::IsValidTaxonMap() const {
  std::vector<bool> id_exists(TaxonCount());
  // Get all ids from map.
  for (const auto &[name, taxon_id] : dag_taxa_) {
    std::ignore = name;
    if (taxon_id.value_ > id_exists.size()) {
      return false;
    }
    if (id_exists[taxon_id.value_]) {
      return false;
    } else {
      id_exists[taxon_id.value_] = true;
    }
  }
  // Taxon map should cover the all ids from 0 to taxon_count-1.
  for (size_t i = 0; i < id_exists.size(); i++) {
    if (!id_exists[i]) {
      return false;
    }
  }
  return true;
}

// ** Reindexer methods:

// NOTE: To be performed *after* DAG modification.
Reindexer SubsplitDAG::BuildNodeReindexer(const size_t prev_node_count) {
  Reindexer node_reindexer = Reindexer::IdentityReindexer(NodeCount());
  // Begin reindex values at taxon count to account for ...leaves?
  size_t running_traversal_idx = taxon_count_;
  NodeId dag_root_node_id = NodeId(prev_node_count - 1);
  // Build node_reindexer by using post-order traversal (topological sort) of entire
  // DAG to assign new ids, where the index is the "before" node_id (stored in the
  // node object), and the value is the "after" node_id.
  DepthFirstWithAction(
      {dag_root_node_id},
      SubsplitDAGTraversalAction(
          // BeforeNode
          [](NodeId node_id) {},
          // AfterNode
          [&node_reindexer, &running_traversal_idx](NodeId node_id) {
            node_reindexer.SetReindex(node_id.value_, running_traversal_idx);
            running_traversal_idx++;
          },
          // BeforeNodeClade
          [](NodeId node_id, bool is_edge_on_left) {},
          // VisitEdge
          [](NodeId node_id, NodeId child_id, bool is_edge_on_left) {}));
  return node_reindexer;
}

Reindexer SubsplitDAG::BuildEdgeReindexer(const size_t prev_edge_count) {
  Reindexer edge_reindexer = Reindexer::IdentityReindexer(EdgeCountWithLeafSubsplits());
  // Only edges from an existing parent node to a new child node need to be reindexed.
  // See SubsplitDAG::AddNodePair().
  for (EdgeId edge_idx = EdgeId(prev_edge_count);
       edge_idx < EdgeCountWithLeafSubsplits(); edge_idx++) {
    // Find edge with given idx.
    auto element = storage_.GetLine(edge_idx);
    Assert(element.has_value(),
           "An edge with given edge_idx did not exist in "
           "SubsplitDAG::BuildEdgeReindexer.");
    auto [node_pair, idx] = element.value();
    std::ignore = idx;
    const auto &[parent_id, child_id] = std::pair<NodeId, NodeId>(node_pair);
    const Bitset parent_subsplit = GetDAGNode(NodeId(parent_id)).GetBitset();
    const Bitset child_subsplit = GetDAGNode(NodeId(child_id)).GetBitset();
    const auto idx_range = GetChildEdgeRange(
        parent_subsplit, child_subsplit.SubsplitIsLeftChildOf(parent_subsplit));
    // New edge is added to the end of the range.
    const auto new_idx =
        EdgeId(edge_reindexer.GetNewIndexByOldIndex(idx_range.second.value_));
    edge_reindexer.ReassignAndShift(edge_idx.value_, new_idx.value_);
  }
  return edge_reindexer;
}

void SubsplitDAG::RemapNodeIds(const Reindexer &node_reindexer) {
  // no need to reindex if no changes were made
  if (node_reindexer == Reindexer::IdentityReindexer(node_reindexer.size())) {
    return;
  }
  std::vector<DAGVertex> nodes = {storage_.GetVertices().begin(),
                                  storage_.GetVertices().end()};
  std::vector<DAGVertex> nodes_copy = Reindexer::Reindex(nodes, node_reindexer);
  storage_.SetVertices(nodes_copy);

  // Update each node's id and leafward/rootward ids.
  for (NodeId node_id = NodeId(0); node_id < NodeCount(); node_id++) {
    GetDAGNode(node_id).RemapNodeIds(node_reindexer);
  }
  // Update `subsplit_to_id_`.
  for (const auto &[subsplit, node_id] : subsplit_to_id_) {
    subsplit_to_id_.at(subsplit) =
        NodeId(node_reindexer.GetNewIndexByOldIndex(node_id.value_));
  }
  // Update edges.
  for (auto i : storage_.GetLines()) {
    storage_.ReindexLine(
        i.GetId(), NodeId(node_reindexer.GetNewIndexByOldIndex(i.GetParent().value_)),
        NodeId(node_reindexer.GetNewIndexByOldIndex(i.GetChild().value_)));
  }
}

void SubsplitDAG::RemapEdgeIdxs(const Reindexer &edge_reindexer) {
  // no need to reindex if no changes were made
  if (edge_reindexer == Reindexer::IdentityReindexer(edge_reindexer.size())) {
    return;
  }
  // Update edges.
  std::vector<DAGLineStorage> edges_copy(storage_.GetLines().size());
  for (auto i : storage_.GetLines()) {
    EdgeId new_idx = EdgeId(edge_reindexer.GetNewIndexByOldIndex(i.GetId().value_));
    edges_copy[new_idx.value_] = i;
    edges_copy[new_idx.value_].SetId(new_idx);
  }
  storage_.SetLines(edges_copy);
  for (NodeId node_id = NodeId(0); node_id < NodeCount(); node_id++) {
    GetDAGNode(node_id).RemapEdgeIdxs(edge_reindexer);
  }
  // Update `parent_to_child_range_`.
  for (const auto &[subsplit, idx_range] : parent_to_child_range_) {
    parent_to_child_range_.at(subsplit) = {
        EdgeId(edge_reindexer.GetNewIndexByOldIndex(idx_range.first.value_)),
        EdgeId(edge_reindexer.GetNewIndexByOldIndex(idx_range.second.value_))};
  }
}
