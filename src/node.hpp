// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// The Node class is how we express tree topologies.
//
// Nodes are immutable after construction except for the id_ and the leaves_.
// The id_ is provided for applications where it is useful to have the edges
// numbered with a contiguous set of integers. The leaves get
// their indices (which are contiguously numbered from 0 through the leaf
// count minus 1) and the rest get ordered according to a postorder traversal.
// Thus the root always has id equal to the number of nodes in the tree.
// See ExampleTopologies below for some examples.
//
// Because this integer assignment cannot be known as we
// are building up the tree, we must make a second pass through the tree, which
// must mutate state. However, this re-id-ing pass is itself deterministic, so
// doing it a second time will always give the same result.
//
// leaves_ is a bitset indicating the set of leaves below. Similarly it needs to
// be calculated on a second pass, because we don't even know the size of the
// bitset as the tree is being built.
//
// Both of these features are prepared using the Polish method.
//
// In summary, call Polish after building your tree if you need to use internal
// node ids or leaf sets. Note that Tree construction calls Polish, if you are
// manually manipulating the topology make you do manipulations with that in
// mind.
//
// Equality is in terms of tree topologies. These mutable members don't matter.

#pragma once

#include <functional>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "bitset.hpp"
#include "intpack.hpp"
#include "sugar.hpp"

class Node {
 public:
  using NodePtr = std::shared_ptr<Node>;
  using Topology = NodePtr;
  using NodePtrVec = std::vector<NodePtr>;
  using NodePtrVecPtr = std::shared_ptr<NodePtrVec>;
  using TopologyCounter = std::unordered_map<NodePtr, uint32_t>;
  // This is the type of functions that are used in the PCSP recursion
  // functions. See `doc/svg/pcsp.svg` for a diagram of the PCSP traversal. In that
  // file, the first tree shows the terminology, and the subsequent trees show
  // the calls to f_root and f_internal.
  //
  // The signature is in 5 parts. The first 4 describe the position in the tree
  // and then the direction: the sister clade, the focal clade, child 0, and
  // child 1. False means down the tree structure and true means up. The 5th
  // part is the top of the virtual root clade, namely the clade containing the
  // virtual root (shown in gray in the diagram). Caution: in the case where the
  // virtual root clade is above the subsplit, the "virtual root clade" will be
  // the entire tree. There's nothing else we can do without rerooting the tree.
  // It's not too hard to exclude the undesired bits with a conditional tree
  // traversal. See IndexerRepresentationOfTopology for an example.
  using UnrootedPCSPFun =
      std::function<void(const Node*, bool, const Node*, bool, const Node*, bool,
                         const Node*, bool, const Node*)>;
  // The rooted version just uses: sister clade, the focal clade, child 0, and child 1.
  using RootedPCSPFun =
      std::function<void(const Node*, const Node*, const Node*, const Node*)>;
  using TwoNodeFun = std::function<void(const Node*, const Node*)>;
  // A function that takes the following node arguments: grandparent, parent, sister,
  // child0, child1.
  using NeighborFun = std::function<void(const Node*, const Node*, const Node*,
                                         const Node*, const Node*)>;

 public:
  explicit Node(uint32_t leaf_id, Bitset leaves);
  explicit Node(NodePtrVec children, size_t id, Bitset leaves);

  size_t Id() const { return id_; }
  uint64_t Tag() const { return tag_; }
  const Bitset& Leaves() const { return leaves_; }
  std::string TagString() const { return StringOfPackedInt(this->tag_); }
  uint32_t MaxLeafID() const { return MaxLeafIDOfTag(tag_); }
  uint32_t LeafCount() const { return LeafCountOfTag(tag_); }
  size_t Hash() const { return hash_; }
  bool IsLeaf() const { return children_.empty(); }
  const NodePtrVec& Children() const { return children_; }

  // Creates a subsplit bitset from given node. Requires tree must be bifurcating.
  Bitset BuildSubsplit() const;
  // Creates an edge PCSP from edge below given clade's side. Requires tree must be
  // bifurcating.
  Bitset BuildPCSP(const SubsplitClade clade) const;

  // Creates a vector of all subsplit bitsets for all nodes in topology.
  std::vector<Bitset> BuildVectorOfSubsplits() const;
  // Creates a vector of all PCSP bitsets for all edges in topology.
  std::vector<Bitset> BuildVectorOfPCSPs() const;

  bool operator==(const Node& other) const;

  NodePtr DeepCopy() const;

  void Preorder(std::function<void(const Node*)> f) const;
  // ConditionalPreorder continues to recur as long as f returns true.
  void ConditionalPreorder(std::function<bool(const Node*)> f) const;
  void Postorder(std::function<void(const Node*)> f) const;
  void LevelOrder(std::function<void(const Node*)> f) const;
  // Apply the pre function before recurring down the tree, and then apply the
  // post function as we are recurring back up the tree.
  void DepthFirst(std::function<void(const Node*)> pre,
                  std::function<void(const Node*)> post) const;

  // We take in two functions, f_root, and f_internal, each of which take three
  // edges.
  // We assume that f_root is symmetric in its last two arguments so that
  // f_root's signature actually looks like f_root(node0, {node1, node2}).
  // We apply f_root to the descendant edges like so: 012, 120, and 201. Because
  // f_root is symmetric in the last two arguments, we are going to get all of
  // the distinct calls of f.
  // At the internal nodes we cycle through triples of (node, sister, parent)
  // for f_internal.
  void TriplePreorder(
      std::function<void(const Node*, const Node*, const Node*)> f_root,
      std::function<void(const Node*, const Node*, const Node*)> f_internal) const;
  // Iterate f through (node, sister, parent) for bifurcating trees using a
  // preorder traversal.
  void TriplePreorderBifurcating(
      std::function<void(const Node*, const Node*, const Node*)> f) const;
  // As above, but getting indices rather than nodes themselves.
  void TripleIdPreorderBifurcating(std::function<void(size_t, size_t, size_t)> f) const;

  // These two functions take functions accepting triples of (node_id,
  // child0_id, child1_id) and apply them according to various traversals.
  void BinaryIdPreorder(std::function<void(size_t, size_t, size_t)> f) const;
  void BinaryIdPostorder(std::function<void(size_t, size_t, size_t)> f) const;

  // See the typedef of UnrootedPCSPFun and RootedPCSPFun to understand the argument
  // type to these functions.
  void UnrootedPCSPPreorder(UnrootedPCSPFun f) const;
  // Apply a RootedPCSPFun to the nodes through a preorder traversal. When allow_leaves
  // is on, the function will be applied on both the internal and leaf nodes.
  // Otherwise, it is only applied on the internal nodes.
  void RootedPCSPPreorder(RootedPCSPFun f, bool allow_leaves) const;
  // Iterate over (leaf sister, leaf) pairs in order. Rooted because that's the only
  // case in which we are guaranteed to have a well defined set of such pairs.
  void RootedSisterAndLeafTraversal(TwoNodeFun f) const;

  // This function prepares the id_ and leaves_ member variables as described at
  // the start of this document. It returns a map that maps the tags to their
  // indices. It's the verb, not the nationality.
  TagSizeMap Polish();

  NodePtr Deroot();

  // ** I/O

  // Return a vector such that the ith component describes the indices for nodes
  // above the current node.
  SizeVectorVector IdsAbove() const;
  // Build a map from each node's id to its parent node. For rootward traversal of tree.
  std::unordered_map<size_t, const Node*> BuildParentNodeMap() const;

  // Output as Newick string, with option for branch lengths.
  std::string Newick(std::function<std::string(const Node*)> node_labeler,
                     const DoubleVectorOption& branch_lengths = std::nullopt) const;
  // Output as Newick string, with options for branch lengths and labels.
  std::string Newick(const DoubleVectorOption& branch_lengths = std::nullopt,
                     const TagStringMapOption& node_labels = std::nullopt,
                     bool show_tags = false) const;

  // Construct a vector such that the ith entry is the id of the parent of the
  // node having id i. We assume that the indices are contiguous, and that the
  // root has the largest id.
  std::vector<size_t> ParentIdVector() const;

  // Outputs this node's id, adjacent leaf ids, and leaf clade bitset to string.
  std::string NodeIdAndLeavesToString() const;
  // Outputs `NodeIdAndLeavesToString` for all entire topology below this node.
  std::string NodeIdAndLeavesToStringForTopology() const;

  // ** Static methods

  // Constructs a leaf node with given id, and an empty taxon clade by default for
  // its leaves.
  static NodePtr Leaf(uint32_t id, Bitset leaves = Bitset(0));
  // Constructs a leaf node with given id, and a single taxon clade with a length of
  // taxon_count for its leaves.
  static NodePtr Leaf(uint32_t id, size_t taxon_count);
  // Join builds a Node with the given descendants, or-ing the leaves_ of the
  // descendants.
  static NodePtr Join(NodePtrVec children, size_t id = SIZE_MAX);
  static NodePtr Join(NodePtr left, NodePtr right, size_t id = SIZE_MAX);
  // Build a tree given a vector of indices, such that each entry gives the
  // id of its parent. We assume that the indices are contiguous, and that
  // the root has the largest id.
  static NodePtr OfParentIdVector(const std::vector<size_t>& indices);

  //     topology           with internal node indices
  //     --------           --------------------------
  // 0: (0,1,(2,3))         (0,1,(2,3)4)5;
  // 1; (0,1,(2,3)) again   (0,1,(2,3)4)5;
  // 2: (0,2,(1,3))         (0,2,(1,3)4)5;
  // 3: (0,(1,(2,3)))       (0,(1,(2,3)4)5)6;
  static NodePtrVec ExampleTopologies();

  // Make a maximally-unbalanced "ladder" tree.
  static NodePtr Ladder(uint32_t leaf_count);

  // A "cryptographic" hash function from Stack Overflow (the std::hash function
  // appears to leave uint32_ts as they are, which doesn't work for our
  // application).
  // https://stackoverflow.com/a/12996028/467327
  static uint32_t SOHash(uint32_t x);

  // Bit rotation from Stack Overflow.
  // c is the amount by which we rotate.
  // https://stackoverflow.com/a/776523/467327
  static size_t SORotate(size_t n, uint32_t c);

 private:
  // Vector of direct child descendants of node in tree topology.
  NodePtrVec children_;
  // NOTE: See beginning of file for notes about the id and the leaves.
  // Unique identifier in tree containing node.
  size_t id_;
  // Bitset of all leaves below node (alternatively can view a leaf as a member of the
  // taxon set in the tree).
  Bitset leaves_;
  // The tag_ is a pair of packed integers representing (1) the maximum leaf ID
  // of the leaves below this node, and (2) the number of leaves below the node.
  uint64_t tag_;
  // Hashkey for node maps.
  size_t hash_;

  // Make copy constructors private to eliminate copying.
  Node(const Node&);
  Node& operator=(const Node&);

  // This is a private Postorder that can change the Node.
  void MutablePostorder(std::function<void(Node*)> f);

  std::string NewickAux(std::function<std::string(const Node*)> node_labeler,
                        const DoubleVectorOption& branch_lengths) const;

  // Make a leaf bitset by or-ing the leaf bitsets of the provided children.
  // Private just to avoid polluting the public interface.
  static Bitset LeavesOf(const NodePtrVec& children);
};

// Compare NodePtrs by their Nodes.
inline bool operator==(const Node::NodePtr& lhs, const Node::NodePtr& rhs) {
  return *lhs == *rhs;
}

inline bool operator!=(const Node::NodePtr& lhs, const Node::NodePtr& rhs) {
  return !(lhs == rhs);
}

namespace std {
template <>
struct hash<Node::NodePtr> {
  size_t operator()(const Node::NodePtr& n) const { return n->Hash(); }
};
template <>
struct equal_to<Node::NodePtr> {
  bool operator()(const Node::NodePtr& lhs, const Node::NodePtr& rhs) const {
    return lhs == rhs;
  }
};
}  // namespace std
