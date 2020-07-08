// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_dag.hpp"
#include "numerical_utils.hpp"

using namespace GPOperations;

void PrintPCSPIndexerFree(const BitsetSizeMap &pcsp_indexer) {
  for (auto it = pcsp_indexer.begin(); it != pcsp_indexer.end(); ++it) {
    std::cout << it->first.SplitChunk(0).ToString() << "|" << it->first.SplitChunk(1).ToString() << ", " << it->second << "\n";
  }
}

GPDAG::GPDAG() : taxon_count_(0), gpcsp_count_(0) {}

GPDAG::GPDAG(const RootedTreeCollection &tree_collection) {
  taxon_count_ = tree_collection.TaxonCount();
  ProcessTrees(tree_collection);
  BuildNodes();
  BuildEdges();
  BuildPCSPIndexer();
  rootward_order_ = RootwardPassTraversal();
  leafward_order_ = LeafwardPassTraversal();
}

size_t GPDAG::NodeCount() const { return dag_nodes_.size(); }
size_t GPDAG::GPCSPCount() const { return gpcsp_count_; }

size_t GPDAG::ContinuousParameterCount() const {
  // Get number of parameters involving fake subsplits.
  size_t fake_subsplit_parameter_count = 0;
  for (size_t i = 0; i < taxon_count_; i++) {
    fake_subsplit_parameter_count += dag_nodes_[i]->GetRootwardRotated().size();
    fake_subsplit_parameter_count += dag_nodes_[i]->GetRootwardSorted().size();
  }

  return GPCSPCount() + fake_subsplit_parameter_count;
}

void GPDAG::ProcessTrees(const RootedTreeCollection &tree_collection) {
  size_t index = 0;
  auto topology_counter = tree_collection.TopologyCounter();
  // Start by adding the rootsplits.
  for (const auto &iter : RootedSBNMaps::RootsplitCounterOf(topology_counter)) {
    SafeInsert(indexer_, iter.first, index);
    rootsplits_.push_back(iter.first);
    index++;
  }
  // Now add the PCSSs.
  for (const auto &[parent, child_counter] :
       RootedSBNMaps::PCSSCounterOf(topology_counter)) {
    SafeInsert(parent_to_range_, parent, {index, index + child_counter.size()});
    for (const auto &child_iter : child_counter) {
      const auto &child = child_iter.first;
      SafeInsert(indexer_, parent + child, index);
      SafeInsert(index_to_child_, index, Bitset::ChildSubsplit(parent, child));
      index++;
    }
  }
  gpcsp_count_ = index;
}

void GPDAG::CreateAndInsertNode(const Bitset &subsplit) {
  if (!subsplit_to_index_.count(subsplit)) {
    size_t id = dag_nodes_.size();
    subsplit_to_index_[subsplit] = id;
    dag_nodes_.push_back(std::make_shared<GPDAGNode>(id, subsplit));
  }
}

void GPDAG::ConnectNodes(size_t idx, bool rotated) {
  auto node = dag_nodes_[idx];
  // Retrieve children subsplits, set edge relation.
  Bitset subsplit = rotated ? node->GetBitset().RotateSubsplit() : node->GetBitset();
  auto children = GetChildrenSubsplits(subsplit, true);
  for (auto child_subsplit : children) {
    auto child_node = dag_nodes_[subsplit_to_index_[child_subsplit]];
    if (rotated) {
      node->AddLeafwardRotated(child_node->Id());
      child_node->AddRootwardRotated(node->Id());
    } else {
      node->AddLeafwardSorted(child_node->Id());
      child_node->AddRootwardSorted(node->Id());
    }
  }
}

// This function returns empty vector if subsplit is invalid or has no child.
std::vector<Bitset> GPDAG::GetChildrenSubsplits(const Bitset &subsplit,
                                                bool include_fake_subsplits) {
  std::vector<Bitset> children_subsplits;

  if (parent_to_range_.count(subsplit)) {
    auto range = parent_to_range_.at(subsplit);
    for (auto idx = range.first; idx < range.second; idx++) {
      auto child_subsplit = index_to_child_.at(idx);
      children_subsplits.push_back(child_subsplit);
    }
  } else if (include_fake_subsplits) {
    // In the case where second chunk of the subsplit is a trivial subsplit,
    // it will not map to any value (parent_to_range_[subsplit] doesn't exist).
    // But we still need to create and connect to fake subsplits in the DAG.
    // A subsplit has a fake subsplit as a child if the first chunk is non-zero
    // and the second chunk has exactly one bit set to 1.
    if (subsplit.SplitChunk(0).Any() &&
              subsplit.SplitChunk(1).SingletonOption().has_value()) {
      // The fake subsplit corresponds to the second chunk of subsplit.
      // Prepend it by 0's.
      Bitset zero(subsplit.size() / 2);
      Bitset fake_subsplit = zero + subsplit.SplitChunk(1);
      children_subsplits.push_back(fake_subsplit);
    }
  }
  
  return children_subsplits;
}

void GPDAG::BuildNodesDepthFirst(const Bitset &subsplit,
                                 std::unordered_set<Bitset> &visited_subsplits) {
  if (!visited_subsplits.count(subsplit)) {
    visited_subsplits.insert(subsplit);
    auto children_subsplits = GetChildrenSubsplits(subsplit, false);
    for (auto child_subsplit : children_subsplits) {
      BuildNodesDepthFirst(child_subsplit, visited_subsplits);
    }
    children_subsplits = GetChildrenSubsplits(subsplit.RotateSubsplit(),
                                              false);
    for (auto child_subsplit : children_subsplits) {
      BuildNodesDepthFirst(child_subsplit, visited_subsplits);
    }

    CreateAndInsertNode(subsplit);
  }
}

void GPDAG::BuildNodes() {
  std::unordered_set<Bitset> visited_subsplits;

  // We will create fake subsplits and insert to dag_nodes_.
  // These nodes will take IDs in [0, taxon_count_).
  Bitset zero(taxon_count_);
  for (size_t i = 0; i < taxon_count_; i++) {
    Bitset fake(taxon_count_);
    fake.set(i);
    Bitset fake_subsplit = zero + fake;
    CreateAndInsertNode(fake_subsplit);
  }
  // We are going to add the remaining nodes.
  // The root splits will take on the higher IDs compared to the non-rootsplits.
  for (auto rootsplit : rootsplits_) {
    auto subsplit = rootsplit + ~rootsplit;
    BuildNodesDepthFirst(subsplit,
                         visited_subsplits);
  }
}

void GPDAG::BuildEdges() {
  for (size_t i = taxon_count_; i < dag_nodes_.size(); i++) {
    ConnectNodes(i, false);
    ConnectNodes(i, true);
  }
}

void GPDAG::Print() {
  for (size_t i = 0; i < dag_nodes_.size(); i++) {
    std::cout << dag_nodes_[i]->ToString() << std::endl;
  }
}

void GPDAG::PrintPCSPIndexer() { PrintPCSPIndexerFree(pcsp_indexer_); }

EigenVectorXd GPDAG::BuildUniformQ() {
  EigenVectorXd q = EigenVectorXd::Ones(ContinuousParameterCount());
  q.segment(0, rootsplits_.size()).array() = 1. / rootsplits_.size();
  for (auto it = subsplit2range_.begin(); it != subsplit2range_.end(); ++it) {
    auto range = subsplit2range_[it->first];
    auto num_child_subsplits = range.second - range.first;
    double val = 1. / num_child_subsplits;
    q.segment(range.first, num_child_subsplits).array() = val;
  }
  return q;
}

GPOperationVector GPDAG::SetRhatToStationary() {
  // Set rhat(s) = stationary for the rootsplits s.
  std::vector<GPOperation> operations;
  for (size_t i = 0; i < rootsplits_.size(); i++) {
    auto rootsplit = rootsplits_[i];
    auto root_idx = subsplit_to_index_[rootsplit + ~rootsplit];
    operations.push_back(SetToStationaryDistribution{
        GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), root_idx), i});
  }
  return operations;
}

void RootwardDepthFirst(size_t id, std::vector<std::shared_ptr<GPDAGNode>> &dag_nodes,
                        std::vector<size_t> &visit_order,
                        std::unordered_set<size_t> &visited_nodes) {
  visited_nodes.insert(id);
  for (size_t child_id : dag_nodes.at(id)->GetRootwardSorted()) {
    if (!visited_nodes.count(child_id)) {
      RootwardDepthFirst(child_id, dag_nodes, visit_order, visited_nodes);
    }
  }
  for (size_t child_id : dag_nodes.at(id)->GetRootwardRotated()) {
    if (!visited_nodes.count(child_id)) {
      RootwardDepthFirst(child_id, dag_nodes, visit_order, visited_nodes);
    }
  }
  visit_order.push_back(id);
}

void LeafwardDepthFirst(size_t id, std::vector<std::shared_ptr<GPDAGNode>> &dag_nodes,
                        std::vector<size_t> &visit_order,
                        std::unordered_set<size_t> &visited_nodes) {
  visited_nodes.insert(id);
  for (size_t child_id : dag_nodes.at(id)->GetLeafwardSorted()) {
    if (!visited_nodes.count(child_id)) {
      LeafwardDepthFirst(child_id, dag_nodes, visit_order, visited_nodes);
    }
  }
  for (size_t child_id : dag_nodes.at(id)->GetLeafwardRotated()) {
    if (!visited_nodes.count(child_id)) {
      LeafwardDepthFirst(child_id, dag_nodes, visit_order, visited_nodes);
    }
  }
  visit_order.push_back(id);
}

std::vector<size_t> GPDAG::LeafwardPassTraversal() {
  // TODO Can we re-add this assert?
  //  Assert(leaf_idx < taxon_count_, std::to_string(leaf_idx) + " >= num taxa.");
  std::vector<size_t> visit_order;
  std::unordered_set<size_t> visited_nodes;
  for (size_t leaf_idx = 0; leaf_idx < taxon_count_; leaf_idx++) {
    RootwardDepthFirst(leaf_idx, dag_nodes_, visit_order, visited_nodes);
  }

  return visit_order;
}

std::vector<size_t> GPDAG::RootwardPassTraversal() {
  std::vector<size_t> visit_order;
  std::unordered_set<size_t> visited_nodes;
  for (auto rootsplit : rootsplits_) {
    size_t root_idx = subsplit_to_index_[rootsplit + ~rootsplit];
    LeafwardDepthFirst(root_idx, dag_nodes_, visit_order, visited_nodes);
  }
  return visit_order;
}

void GPDAG::BuildPCSPIndexer() {
  size_t idx = 0;
  for (auto rootsplit : rootsplits_) {
    SafeInsert(pcsp_indexer_, rootsplit + ~rootsplit, idx);
    idx++;
  }

  for (size_t i = taxon_count_; i < dag_nodes_.size(); i++) {
    auto node = dag_nodes_[i];
    auto n_child = node->GetLeafwardSorted().size();
    if (n_child > 0) {
      SafeInsert(subsplit2range_, node->GetBitset(), {idx, idx + n_child});
      for (size_t j = 0; j < n_child; j++) {
        auto child = dag_nodes_[node->GetLeafwardSorted()[j]];
        SafeInsert(pcsp_indexer_, node->GetBitset() + child->GetBitset(), idx);
        idx++;
      }
    }
    n_child = node->GetLeafwardRotated().size();
    if (n_child > 0) {
      SafeInsert(subsplit2range_, node->GetBitset().RotateSubsplit(),
                 {idx, idx + n_child});
      for (size_t j = 0; j < n_child; j++) {
        auto child = dag_nodes_[node->GetLeafwardRotated()[j]];
        SafeInsert(pcsp_indexer_,
                   node->GetBitset().RotateSubsplit() + child->GetBitset(), idx);
        idx++;
      }
    }
  }
}

size_t GetPLVIndex(PLVType plv_type, size_t node_count, size_t src_idx) {
  switch (plv_type) {
    case PLVType::P:
      return src_idx;
    case PLVType::P_HAT:
      return node_count + src_idx;
    case PLVType::P_HAT_TILDE:
      return 2 * node_count + src_idx;
    case PLVType::R_HAT:
      return 3 * node_count + src_idx;
    case PLVType::R:
      return 4 * node_count + src_idx;
    case PLVType::R_TILDE:
      return 5 * node_count + src_idx;
    default:
      Failwith("Invalid PLV index requested.");
  }
}

void GPDAG::AddRootwardWeightedSumAccumulateOperations(std::shared_ptr<GPDAGNode> node,
                                                       bool rotated,
                                                       GPOperationVector &operations) {
  std::vector<size_t> child_idxs =
      rotated ? node->GetLeafwardRotated() : node->GetLeafwardSorted();
  PLVType plv_type = rotated ? PLVType::P_HAT_TILDE : PLVType::P_HAT;
  auto parent_subsplit =
      rotated ? node->GetBitset().RotateSubsplit() : node->GetBitset();
  for (size_t child_idx : child_idxs) {
    auto child_node = dag_nodes_[child_idx];
    auto child_subsplit = child_node->GetBitset();
    auto pcsp = parent_subsplit + child_subsplit;
    if (!pcsp_indexer_.count(pcsp)) {
      Failwith("Non-existent PCSP index.");
    }
    auto pcsp_idx = pcsp_indexer_[pcsp];

    operations.push_back(WeightedSumAccumulate{
        GetPLVIndex(plv_type, dag_nodes_.size(), node->Id()), pcsp_idx,
        GetPLVIndex(PLVType::P, dag_nodes_.size(), child_idx)});
  }
}

void GPDAG::AddLeafwardWeightedSumAccumulateOperations(std::shared_ptr<GPDAGNode> node,
                                                       GPOperationVector &operations) {
  auto subsplit = node->GetBitset();
  for (size_t parent_idx : node->GetRootwardSorted()) {
    auto parent_node = dag_nodes_[parent_idx];
    auto parent_subsplit = parent_node->GetBitset();
    auto pcsp_idx = pcsp_indexer_[parent_subsplit + subsplit];

    operations.push_back(WeightedSumAccumulate{
        GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), node->Id()), pcsp_idx,
        GetPLVIndex(PLVType::R, dag_nodes_.size(), parent_node->Id())});
  }
  for (size_t parent_idx : node->GetRootwardRotated()) {
    auto parent_node = dag_nodes_[parent_idx];
    auto parent_subsplit = parent_node->GetBitset().RotateSubsplit();
    auto pcsp_idx = pcsp_indexer_[parent_subsplit + subsplit];

    operations.push_back(WeightedSumAccumulate{
        GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), node->Id()), pcsp_idx,
        GetPLVIndex(PLVType::R_TILDE, dag_nodes_.size(), parent_node->Id())});
  }
}

void GPDAG::OptimizeSBNParameters(const Bitset &subsplit,
                                  GPOperationVector &operations) {
  if (subsplit2range_.count(subsplit)) {
    auto param_range = subsplit2range_[subsplit];
    if (param_range.second - param_range.first > 1) {
      operations.push_back(
          UpdateSBNProbabilities{param_range.first, param_range.second});
    }
  }
}

size_t GPDAG::GetPCSPIndex(size_t parent_node_idx, size_t child_node_idx,
                           bool rotated) {
  auto parent_node = dag_nodes_[parent_node_idx];
  auto child_node = dag_nodes_[child_node_idx];
  auto pcsp = rotated
                  ? parent_node->GetBitset().RotateSubsplit() + child_node->GetBitset()
                  : parent_node->GetBitset() + child_node->GetBitset();
  return pcsp_indexer_[pcsp];
}

GPOperationVector GPDAG::RootwardPass(std::vector<size_t> visit_order) {
  // Perform first rootward pass. No optimization.
  GPOperationVector operations;
  for (size_t node_idx : visit_order) {
    auto node = dag_nodes_[node_idx];
    if (node->IsLeaf()) {
      continue;
    }

    // Perform the following operations:
    // 1. WeightedSumAccummulate to get phat(s), phat(s_tilde).
    // 2. Multiply to get p(s) = phat(s) \circ phat(s_tilde).
    AddRootwardWeightedSumAccumulateOperations(node, false, operations);
    AddRootwardWeightedSumAccumulateOperations(node, true, operations);
    operations.push_back(
        Multiply{node_idx, GetPLVIndex(PLVType::P_HAT, dag_nodes_.size(), node_idx),
                 GetPLVIndex(PLVType::P_HAT_TILDE, dag_nodes_.size(), node_idx)});
  }

  return operations;
}

GPOperationVector GPDAG::RootwardPass() { return RootwardPass(rootward_order_); }

GPOperationVector GPDAG::LeafwardPass(std::vector<size_t> visit_order) {
  GPOperationVector operations;
  for (size_t node_idx : visit_order) {
    auto node = dag_nodes_[node_idx];

    // Perform the following operations:
    // 1. WeightedSumAccumulate: rhat(s) += \sum_t q(s|t) P'(s|t) r(t)
    // 2. Multiply: r(s) = rhat(s) \circ phat(s_tilde).
    // 3. Multiply: r(s_tilde) = rhat(s) \circ phat(s).
    AddLeafwardWeightedSumAccumulateOperations(node, operations);
    operations.push_back(
        Multiply{GetPLVIndex(PLVType::R, dag_nodes_.size(), node_idx),
                 GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), node_idx),
                 GetPLVIndex(PLVType::P_HAT_TILDE, dag_nodes_.size(), node_idx)});
    operations.push_back(
        Multiply{GetPLVIndex(PLVType::R_TILDE, dag_nodes_.size(), node_idx),
                 GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), node_idx),
                 GetPLVIndex(PLVType::P_HAT, dag_nodes_.size(), node_idx)});
  }
  return operations;
}

GPOperationVector GPDAG::LeafwardPass() { return LeafwardPass(leafward_order_); }

GPOperationVector GPDAG::SetRootwardZero() {
  GPOperationVector operations;
  auto node_count = dag_nodes_.size();
  for (size_t i = taxon_count_; i < node_count; i++) {
    operations.push_back(Zero{GetPLVIndex(PLVType::P, dag_nodes_.size(), i)});
    operations.push_back(Zero{GetPLVIndex(PLVType::P_HAT, dag_nodes_.size(), i)});
    operations.push_back(Zero{GetPLVIndex(PLVType::P_HAT_TILDE, dag_nodes_.size(), i)});
  }
  return operations;
}

GPOperationVector GPDAG::SetLeafwardZero() {
  GPOperationVector operations;
  auto node_count = dag_nodes_.size();
  for (size_t i = 0; i < node_count; i++) {
    operations.push_back(Zero{GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), i)});
    operations.push_back(Zero{GetPLVIndex(PLVType::R, dag_nodes_.size(), i)});
    operations.push_back(Zero{GetPLVIndex(PLVType::R_TILDE, dag_nodes_.size(), i)});
  }
  for (auto rootsplit : rootsplits_) {
    size_t root_idx = subsplit_to_index_[rootsplit + ~rootsplit];
    operations.push_back(SetToStationaryDistribution{
        GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), root_idx)});
  }
  return operations;
}

void GPDAG::ScheduleBranchLengthOptimization(size_t node_id,
                                             std::unordered_set<size_t> &visited_nodes,
                                             GPOperationVector &operations) {
  visited_nodes.insert(node_id);
  auto node = dag_nodes_[node_id];

  if (!node->IsRoot()) {
    // Compute R_HAT(s) = \sum_{t : s < t} q(s|t) P(s|t) r(t).
    operations.push_back(Zero{GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), node_id)});
    for (size_t parent_id : node->GetRootwardSorted()) {
      auto parent_node = dag_nodes_[parent_id];
      size_t pcsp_idx = pcsp_indexer_[parent_node->GetBitset() + node->GetBitset()];
      operations.push_back(WeightedSumAccumulate{
          GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), node_id), pcsp_idx,
          GetPLVIndex(PLVType::R, dag_nodes_.size(), parent_id)});
    }
    for (size_t parent_id : node->GetRootwardRotated()) {
      auto parent_node = dag_nodes_[parent_id];
      size_t pcsp_idx =
          pcsp_indexer_[parent_node->GetBitset().RotateSubsplit() + node->GetBitset()];
      operations.push_back(WeightedSumAccumulate{
          GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), node_id), pcsp_idx,
          GetPLVIndex(PLVType::R_TILDE, dag_nodes_.size(), parent_id)});
    }

    // Update r(s) and r_tilde(s)
    operations.push_back(
        Multiply{GetPLVIndex(PLVType::R, dag_nodes_.size(), node_id),
                 GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), node_id),
                 GetPLVIndex(PLVType::P_HAT_TILDE, dag_nodes_.size(), node_id)});
    operations.push_back(
        Multiply{GetPLVIndex(PLVType::R_TILDE, dag_nodes_.size(), node_id),
                 GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), node_id),
                 GetPLVIndex(PLVType::P_HAT, dag_nodes_.size(), node_id)});
  }

  if (node->IsLeaf()) return;

  operations.push_back(Zero{GetPLVIndex(PLVType::P_HAT, dag_nodes_.size(), node_id)});
  for (size_t child_id : dag_nodes_.at(node_id)->GetLeafwardSorted()) {
    if (!visited_nodes.count(child_id)) {
      ScheduleBranchLengthOptimization(child_id, visited_nodes, operations);
    }

    auto child_node = dag_nodes_[child_id];
    size_t pcsp_idx = pcsp_indexer_[node->GetBitset() + child_node->GetBitset()];
    operations.push_back(OptimizeBranchLength{
        GetPLVIndex(PLVType::P, dag_nodes_.size(), child_id),
        GetPLVIndex(PLVType::R, dag_nodes_.size(), node_id), pcsp_idx});
    // Update p_hat(s)
    operations.push_back(WeightedSumAccumulate{
        GetPLVIndex(PLVType::P_HAT, dag_nodes_.size(), node_id),
        pcsp_idx,
        GetPLVIndex(PLVType::P, dag_nodes_.size(), child_id),
    });
  }
  // Update r_tilde(t) = r_hat(t) \circ p_hat(t).
  operations.push_back(
      Multiply{GetPLVIndex(PLVType::R_TILDE, dag_nodes_.size(), node_id),
               GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), node_id),
               GetPLVIndex(PLVType::P_HAT, dag_nodes_.size(), node_id)});

  operations.push_back(
      Zero{GetPLVIndex(PLVType::P_HAT_TILDE, dag_nodes_.size(), node_id)});
  for (size_t child_id : dag_nodes_.at(node_id)->GetLeafwardRotated()) {
    if (!visited_nodes.count(child_id)) {
      ScheduleBranchLengthOptimization(child_id, visited_nodes, operations);
    }
    auto child_node = dag_nodes_[child_id];
    size_t pcsp_idx =
        pcsp_indexer_[node->GetBitset().RotateSubsplit() + child_node->GetBitset()];
    operations.push_back(OptimizeBranchLength{
        GetPLVIndex(PLVType::P, dag_nodes_.size(), child_id),
        GetPLVIndex(PLVType::R_TILDE, dag_nodes_.size(), node_id), pcsp_idx});
    operations.push_back(WeightedSumAccumulate{
        GetPLVIndex(PLVType::P_HAT_TILDE, dag_nodes_.size(), node_id),
        pcsp_idx,
        GetPLVIndex(PLVType::P, dag_nodes_.size(), child_id),
    });
  }
  // Update r(t) = r_hat(t) \circ p_hat_tilde(t).
  operations.push_back(
      Multiply{GetPLVIndex(PLVType::R, dag_nodes_.size(), node_id),
               GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), node_id),
               GetPLVIndex(PLVType::P_HAT_TILDE, dag_nodes_.size(), node_id)});

  // Update p(t).
  operations.push_back(
      Multiply{GetPLVIndex(PLVType::P, dag_nodes_.size(), node_id),
               GetPLVIndex(PLVType::P_HAT, dag_nodes_.size(), node_id),
               GetPLVIndex(PLVType::P_HAT_TILDE, dag_nodes_.size(), node_id)});
};

GPOperationVector GPDAG::BranchLengthOptimization() {
  GPOperationVector operations;
  std::unordered_set<size_t> visited_nodes;
  for (auto rootsplit : rootsplits_) {
    ScheduleBranchLengthOptimization(subsplit_to_index_[rootsplit + ~rootsplit],
                                     visited_nodes, operations);
  }
  return operations;
}

void GPDAG::ScheduleSBNParametersOptimization(size_t node_id,
                                              std::unordered_set<size_t> &visited_nodes,
                                              GPOperationVector &operations) {
  visited_nodes.insert(node_id);
  auto node = dag_nodes_[node_id];

  if (!node->IsRoot()) {
    // We compute R_HAT(s) = \sum_{t : s < t} q(s|t) P(s|t) r(t).
    // This is necessary to reflect changes to r(t) as well as new values
    // for q(s|t).
    // TODO: Code refactoring -- similar code is repeated twice.
    operations.push_back(Zero{GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), node_id)});
    for (size_t parent_id : node->GetRootwardSorted()) {
      auto parent_node = dag_nodes_[parent_id];
      size_t pcsp_idx = pcsp_indexer_[parent_node->GetBitset() + node->GetBitset()];
      operations.push_back(WeightedSumAccumulate{
          GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), node_id), pcsp_idx,
          GetPLVIndex(PLVType::R, dag_nodes_.size(), parent_id)});
    }
    for (size_t parent_id : node->GetRootwardRotated()) {
      auto parent_node = dag_nodes_[parent_id];
      size_t pcsp_idx =
          pcsp_indexer_[parent_node->GetBitset().RotateSubsplit() + node->GetBitset()];
      operations.push_back(WeightedSumAccumulate{
          GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), node_id), pcsp_idx,
          GetPLVIndex(PLVType::R_TILDE, dag_nodes_.size(), parent_id)});
    }

    // Update r(s) and r_tilde(s) using rhat(s).
    operations.push_back(
        Multiply{GetPLVIndex(PLVType::R, dag_nodes_.size(), node_id),
                 GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), node_id),
                 GetPLVIndex(PLVType::P_HAT_TILDE, dag_nodes_.size(), node_id)});
    operations.push_back(
        Multiply{GetPLVIndex(PLVType::R_TILDE, dag_nodes_.size(), node_id),
                 GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), node_id),
                 GetPLVIndex(PLVType::P_HAT, dag_nodes_.size(), node_id)});
  }

  if (node->IsLeaf()) return;

  operations.push_back(Zero{GetPLVIndex(PLVType::P_HAT, dag_nodes_.size(), node_id)});
  for (size_t child_id : dag_nodes_.at(node_id)->GetLeafwardSorted()) {
    if (!visited_nodes.count(child_id)) {
      ScheduleSBNParametersOptimization(child_id, visited_nodes, operations);
    }

    auto child_node = dag_nodes_[child_id];
    size_t pcsp_idx = pcsp_indexer_[node->GetBitset() + child_node->GetBitset()];
    // Update p_hat(s)
    operations.push_back(WeightedSumAccumulate{
        GetPLVIndex(PLVType::P_HAT, dag_nodes_.size(), node_id),
        pcsp_idx,
        GetPLVIndex(PLVType::P, dag_nodes_.size(), child_id),
    });
    operations.push_back(
        Likelihood{pcsp_idx, GetPLVIndex(PLVType::R, dag_nodes_.size(), node->Id()),
                   GetPLVIndex(PLVType::P, dag_nodes_.size(), child_node->Id())});
  }
  OptimizeSBNParameters(node->GetBitset(), operations);

  // Update r_tilde(t) = r_hat(t) \circ p_hat(t).
  operations.push_back(
      Multiply{GetPLVIndex(PLVType::R_TILDE, dag_nodes_.size(), node_id),
               GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), node_id),
               GetPLVIndex(PLVType::P_HAT, dag_nodes_.size(), node_id)});

  operations.push_back(
      Zero{GetPLVIndex(PLVType::P_HAT_TILDE, dag_nodes_.size(), node_id)});
  for (size_t child_id : dag_nodes_.at(node_id)->GetLeafwardRotated()) {
    if (!visited_nodes.count(child_id)) {
      ScheduleSBNParametersOptimization(child_id, visited_nodes, operations);
    }
    auto child_node = dag_nodes_[child_id];
    size_t pcsp_idx =
        pcsp_indexer_[node->GetBitset().RotateSubsplit() + child_node->GetBitset()];
    operations.push_back(WeightedSumAccumulate{
        GetPLVIndex(PLVType::P_HAT_TILDE, dag_nodes_.size(), node_id),
        pcsp_idx,
        GetPLVIndex(PLVType::P, dag_nodes_.size(), child_id),
    });
    operations.push_back(Likelihood{
        pcsp_idx, GetPLVIndex(PLVType::R_TILDE, dag_nodes_.size(), node->Id()),
        GetPLVIndex(PLVType::P, dag_nodes_.size(), child_node->Id())});
  }
  OptimizeSBNParameters(node->GetBitset().RotateSubsplit(), operations);

  // Update r(t) = r_hat(t) \circ p_hat_tilde(t).
  operations.push_back(
      Multiply{GetPLVIndex(PLVType::R, dag_nodes_.size(), node_id),
               GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), node_id),
               GetPLVIndex(PLVType::P_HAT_TILDE, dag_nodes_.size(), node_id)});

  // Update p(t).
  operations.push_back(
      Multiply{GetPLVIndex(PLVType::P, dag_nodes_.size(), node_id),
               GetPLVIndex(PLVType::P_HAT, dag_nodes_.size(), node_id),
               GetPLVIndex(PLVType::P_HAT_TILDE, dag_nodes_.size(), node_id)});
};

GPOperationVector GPDAG::SBNParameterOptimization() {
  GPOperationVector operations;
  std::unordered_set<size_t> visited_nodes;
  for (size_t i = 0; i < rootsplits_.size(); i++) {
    auto rootsplit = rootsplits_[i];
    ScheduleSBNParametersOptimization(subsplit_to_index_[rootsplit + ~rootsplit],
                                      visited_nodes, operations);
    auto node_id = subsplit_to_index_[rootsplit + ~rootsplit];
    operations.push_back(
        MarginalLikelihood{GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), node_id), i,
                           GetPLVIndex(PLVType::P, dag_nodes_.size(), node_id)});
  }
  // Optimize SBN parameters for the rootsplits: at this point, p-vectors are
  // already updated in the call to ScheduleSBNParametersOptimization.
  operations.push_back(UpdateSBNProbabilities{0, rootsplits_.size()});
  return operations;
}

GPOperationVector GPDAG::ComputeLikelihoods() {
  // Compute likelihood values l(s|t) for each PCSP s|t.
  GPOperationVector operations;
  for (size_t i = taxon_count_; i < dag_nodes_.size(); i++) {
    auto node = dag_nodes_[i];
    for (size_t child_idx : node->GetLeafwardSorted()) {
      auto child_node = dag_nodes_[child_idx];
      auto pcsp_idx = pcsp_indexer_[node->GetBitset() + child_node->GetBitset()];
      operations.push_back(
          Likelihood{pcsp_idx, GetPLVIndex(PLVType::R, dag_nodes_.size(), node->Id()),
                     GetPLVIndex(PLVType::P, dag_nodes_.size(), child_node->Id())});
    }
    for (size_t child_idx : node->GetLeafwardRotated()) {
      auto child_node = dag_nodes_[child_idx];
      auto pcsp_idx =
          pcsp_indexer_[node->GetBitset().RotateSubsplit() + child_node->GetBitset()];
      operations.push_back(Likelihood{
          pcsp_idx, GetPLVIndex(PLVType::R_TILDE, dag_nodes_.size(), node->Id()),
          GetPLVIndex(PLVType::P, dag_nodes_.size(), child_node->Id())});
    }
  }

  // Compute marginal likelihood.
  for (size_t root_idx = 0; root_idx < rootsplits_.size(); root_idx++) {
    auto rootsplit = rootsplits_[root_idx];
    auto node_idx = subsplit_to_index_[rootsplit + ~rootsplit];
    operations.push_back(MarginalLikelihood{
        GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), node_idx), root_idx,
        GetPLVIndex(PLVType::P, dag_nodes_.size(), node_idx)});
  }

  std::cout << operations;

  return operations;
}

GPOperationVector GPDAG::MarginalLikelihoodOperations() {
  // Compute marginal likelihood.
  GPOperationVector operations;
  for (size_t i = 0; i < rootsplits_.size(); i++) {
    auto rootsplit = rootsplits_[i];
    auto root_subsplit = rootsplit + ~rootsplit;
    size_t root_idx = subsplit_to_index_[root_subsplit];
    operations.push_back(
        MarginalLikelihood{GetPLVIndex(PLVType::R_HAT, dag_nodes_.size(), root_idx), i,
                           GetPLVIndex(PLVType::P, dag_nodes_.size(), root_idx)});
  }
  return operations;
}
