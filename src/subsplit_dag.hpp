// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// The purpose of this class is to hold a DAG built from the parent-child relationships
// of the subsplits.

#ifndef SRC_SUBSPLIT_DAG_HPP_
#define SRC_SUBSPLIT_DAG_HPP_

#include "gp_dag_node.hpp"
#include "rooted_tree_collection.hpp"
#include "sbn_maps.hpp"

class SubsplitDAG {
 public:
  // NodeLambda is for iterating over nodes.
  using NodeLambda = std::function<void(const GPDAGNode *)>;
  // EdgeDestinationLambda takes in a rotation status (true is rotated, false is not)
  // and a "destination" node. For iterating over DAG edges with a rotation status.
  using EdgeDestinationLambda = std::function<void(bool, const GPDAGNode *)>;

  SubsplitDAG();
  explicit SubsplitDAG(const RootedTreeCollection &tree_collection);

  size_t NodeCount() const;
  // How many topologies can be expressed by the GPDAG? Expressed as a double because
  // this number can be big.
  double TopologyCount() const;
  // Each node in a topology is constructed with GPDAGNode ID as Node ID.
  Node::NodePtrVec GenerateAllGPNodeIndexedTopologies() const;
  size_t RootsplitCount() const;
  size_t GPCSPCount() const;
  size_t GPCSPCountWithFakeSubsplits() const;

  void Print() const;
  void PrintGPCSPIndexer() const;

  const BitsetSizeMap &GetGPCSPIndexer() const;
  const BitsetSizePairMap &GetParentToRange() const;
  GPDAGNode *GetDagNode(size_t node_id) const;

  // Access the value of the gpcsp_indexer_ at the given "expanded" PCSP.
  // Asserts to make sure that the PCSP is well formed.
  size_t GetGPCSPIndex(const Bitset &parent_subsplit,
                       const Bitset &child_subsplit) const;

  // Discrete uniform distribution over each subsplit.
  [[nodiscard]] EigenVectorXd BuildUniformQ() const;
  // Uniform prior over all topologies.
  [[nodiscard]] EigenVectorXd BuildUniformPrior() const;

 protected:
  size_t taxon_count_;
  size_t gpcsp_count_without_fake_subsplits_;
  // The collection of rootsplits, with the same indexing as in the indexer_.
  BitsetVector rootsplits_;
  // A map going from the index of a PCSP to its child.
  SizeBitsetMap index_to_child_;
  // This indexer is an expanded version of parent_to_range_ in sbn_instance:
  // it includes single element range for fake subsplits.
  BitsetSizePairMap parent_to_range_;

  // This indexer is a similarly expanded version of indexer_ in an SBNSupport, but also
  // encoding rootsplits as full splits rather than a binary indicator.
  //
  // This indexer is used for q_, branch_lengths_, log_likelihoods_ in GPEngine.
  BitsetSizeMap gpcsp_indexer_;

  // We will call the index of DAG nodes "ids" to distinguish them from GPCSP indexes.
  // This corresponds to the analogous concept for topologies.
  //
  // A map from Bitset to the corresponding index in dag_nodes_.
  // The first entries are reserved for fake subsplits.
  // The last entries are reserved for rootsplits.
  BitsetSizeMap subsplit_to_id_;
  std::vector<std::unique_ptr<GPDAGNode>> dag_nodes_;

  // Total number of topologies spanned by the DAG.
  double topology_count_;
  // Storage for the number of topologies below for each node.
  EigenVectorXd topology_count_below_;

  // Iterate over the "real" nodes, i.e. those that do not correspond to fake subsplits.
  void IterateOverRealNodes(const NodeLambda &f) const;
  // Iterate over the all leafward edges, rotated and sorted, of node using an
  // EdgeDestinationLambda.
  void IterateOverLeafwardEdges(const GPDAGNode *node,
                                const EdgeDestinationLambda &f) const;
  // Iterate over only the rotated/sorted leafward edges of node using a NodeLambda.
  void IterateOverLeafwardEdges(const GPDAGNode *node, bool rotated,
                                const NodeLambda &f) const;
  // Iterate over the rootward edges of node using an EdgeDestinationLambda.
  void IterateOverRootwardEdges(const GPDAGNode *node,
                                const EdgeDestinationLambda &f) const;
  // Iterate over the node ids corresponding to rootsplits.
  void IterateOverRootsplitIds(const std::function<void(size_t)> &f) const;

  // Gives the children subsplits of a given parent subsplit, optionally including fake
  // subsplits.
  std::vector<Bitset> GetChildrenSubsplits(const Bitset &subsplit,
                                           bool include_fake_subsplits = false);

  void ProcessTrees(const RootedTreeCollection &tree_collection);
  void CreateAndInsertNode(const Bitset &subsplit);
  // Connect the `idx` node to its children, and its children to it, rotating as needed.
  void ConnectNodes(size_t idx, bool rotated);

  void BuildNodesDepthFirst(const Bitset &subsplit,
                            std::unordered_set<Bitset> &visited_subsplits);
  void BuildNodes();
  void BuildEdges();
  void CountTopologies();
  // Expand gpcsp_indexer_ and parent_to_range_ with fake subsplits at the end.
  void AddFakeSubsplitsToGPCSPIndexerAndParentToRange();
  // Update gpcsp_indexer_ rootsplits to be full subsplits.
  void ExpandRootsplitsInIndexer();

  [[nodiscard]] SizeVector LeafwardPassTraversal() const;
  [[nodiscard]] SizeVector RootwardPassTraversal() const;

  Bitset PerhapsRotateSubsplit(const Bitset &subsplit, bool rotated);
};

#endif  // SRC_SUBSPLIT_DAG_HPP_
