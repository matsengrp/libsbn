// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.

#include "node.hpp"

#include <algorithm>
#include <climits>
#include <deque>
#include <functional>
#include <iostream>
#include <memory>
#include <sstream>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>

// We assume that the sizes of things related to trees are smaller than
// UINT32_MAX and so use that as our fundamental type. Because we use the STL,
// we use size_t as well, the size of which is implementation-dependent. Here we
// make sure that size_t is big enough (I know this is a little silly because
// size_t is quite big on 64 bit systems, but still...).
static_assert(UINT32_MAX <= SIZE_MAX, "size_t is too small.");

Node::Node(uint32_t leaf_id, Bitset leaves)
    : children_({}),
      id_(leaf_id),
      leaves_(std::move(leaves)),
      tag_(PackInts(leaf_id, 1)),
      hash_(SOHash(leaf_id)) {}

Node::Node(NodePtrVec children, size_t id, Bitset leaves)
    : children_(children), id_(id), leaves_(leaves) {
  Assert(!children_.empty(), "Called internal Node constructor with no children.");
  // Order the children by their max leaf ids.
  std::sort(children_.begin(), children_.end(), [](const auto& lhs, const auto& rhs) {
    if (lhs->MaxLeafID() == rhs->MaxLeafID()) {
      // Children should have non-overlapping leaf sets, so there
      // should not be ties.
      Failwith("Tie observed between " + lhs->Newick() + " and " + rhs->Newick() +
               "\n" + "Do you have a taxon name repeated?");
    }
    return lhs->MaxLeafID() < rhs->MaxLeafID();
  });
  // Children are sorted by their max_leaf_id, so we can get the max by
  // looking at the last element.
  uint32_t max_leaf_id = children_.back()->MaxLeafID();
  uint32_t leaf_count = 0;
  hash_ = 0;
  for (const auto& child : children_) {
    leaf_count += child->LeafCount();
    hash_ ^= child->Hash();
  }
  tag_ = PackInts(max_leaf_id, leaf_count);
  // Bit rotation is necessary because if we only XOR then we can get
  // collisions when identical tips are in different
  // ordered subtrees (an example is in below doctest).
  hash_ = SORotate(hash_, 1);
}

bool Node::operator==(const Node& other) const {
  if (this->Hash() != other.Hash()) {
    return false;
  }
  size_t child_count = this->Children().size();
  if (child_count != other.Children().size()) {
    return false;
  }
  for (size_t i = 0; i < child_count; i++) {
    if (!(*children_[i] == *other.Children()[i])) {
      return false;
    }
  }
  return true;
}

Node::NodePtr Node::DeepCopy() const {
  return Node::OfParentIdVector(ParentIdVector());
}

void Node::Preorder(std::function<void(const Node*)> f) const {
  std::stack<const Node*> stack;
  stack.push(this);
  const Node* node;
  while (stack.size()) {
    node = stack.top();
    stack.pop();
    f(node);
    const auto& children = node->Children();
    for (auto iter = children.rbegin(); iter != children.rend(); ++iter) {
      stack.push((*iter).get());
    }
  }
}

void Node::ConditionalPreorder(std::function<bool(const Node*)> f) const {
  std::stack<const Node*> stack;
  stack.push(this);
  const Node* node;
  while (stack.size()) {
    node = stack.top();
    stack.pop();
    if (f(node)) {
      const auto& children = node->Children();
      for (auto iter = children.rbegin(); iter != children.rend(); ++iter) {
        stack.push((*iter).get());
      }
    }
  }
}

void Node::MutablePostorder(std::function<void(Node*)> f) {
  // The stack records the nodes and whether they have been visited or not.
  std::stack<std::pair<Node*, bool>> stack;
  stack.push({this, false});
  Node* node;
  bool visited;
  while (stack.size()) {
    std::tie(node, visited) = stack.top();
    stack.pop();
    if (visited) {
      // If we've already visited this node then we are on our way back.
      f(node);
    } else {
      // If not then we need to push ourself back on the stack (noting that
      // we've been visited)...
      stack.push({node, true});
      // And all of our children, which have not.
      const auto& children = node->Children();
      for (auto iter = children.rbegin(); iter != children.rend(); ++iter) {
        stack.push({(*iter).get(), false});
      }
    }
  }
}

void Node::Postorder(std::function<void(const Node*)> f) const {
  // https://stackoverflow.com/a/56603436/467327
  Node* mutable_this = const_cast<Node*>(this);
  mutable_this->MutablePostorder(f);
}

void Node::DepthFirst(std::function<void(const Node*)> pre,
                      std::function<void(const Node*)> post) const {
  // The stack records the nodes and whether they have been visited or not.
  std::stack<std::pair<const Node*, bool>> stack;
  stack.push({this, false});
  const Node* node;
  bool visited;
  while (stack.size()) {
    std::tie(node, visited) = stack.top();
    stack.pop();
    if (visited) {
      // If we've already visited this node then we are on our way back.
      post(node);
    } else {
      pre(node);
      // If not then we need to push ourself back on the stack (noting that
      // we've been visited)...
      stack.push({node, true});
      // And all of our children, which have not.
      const auto& children = node->Children();
      for (auto iter = children.rbegin(); iter != children.rend(); ++iter) {
        stack.push({(*iter).get(), false});
      }
    }
  }
}

void Node::LevelOrder(std::function<void(const Node*)> f) const {
  std::deque<const Node*> deque = {this};
  while (deque.size()) {
    auto n = deque.front();
    deque.pop_front();
    f(n);
    for (const auto& child : n->children_) {
      deque.push_back(child.get());
    }
  }
}

static std::function<void(const Node*, const Node*, const Node*)> const TripletIdInfix(
    std::function<void(int, int, int)> f) {
  return [&f](const Node* node0, const Node* node1, const Node* node2) {
    f(static_cast<int>(node0->Id()), static_cast<int>(node1->Id()),
      static_cast<int>(node2->Id()));
  };
}

void Node::TripleIdPreorderBifurcating(
    std::function<void(size_t, size_t, size_t)> f) const {
  TriplePreorderBifurcating(TripletIdInfix(f));
}

static std::function<void(const Node*)> const BinaryIdInfix(
    std::function<void(int, int, int)> f) {
  return [&f](const Node* node) {
    if (!node->IsLeaf()) {
      Assert(node->Children().size() == 2, "BinaryIdInfix expects a bifurcating tree.");
      f(static_cast<int>(node->Id()), static_cast<int>(node->Children()[0]->Id()),
        static_cast<int>(node->Children()[1]->Id()));
    }
  };
}

void Node::BinaryIdPreorder(const std::function<void(size_t, size_t, size_t)> f) const {
  Preorder(BinaryIdInfix(f));
}

void Node::BinaryIdPostorder(
    const std::function<void(size_t, size_t, size_t)> f) const {
  Postorder(BinaryIdInfix(f));
}

void Node::TriplePreorder(
    std::function<void(const Node*, const Node*, const Node*)> f_root,
    std::function<void(const Node*, const Node*, const Node*)> f_internal) const {
  Assert(children_.size() == 3,
         "TriplePreorder expects a tree with a trifurcation at the root.");
  f_root(children_[0].get(), children_[1].get(), children_[2].get());
  children_[0]->TriplePreorderBifurcating(f_internal);
  f_root(children_[1].get(), children_[2].get(), children_[0].get());
  children_[1]->TriplePreorderBifurcating(f_internal);
  f_root(children_[2].get(), children_[0].get(), children_[1].get());
  children_[2]->TriplePreorderBifurcating(f_internal);
}

void Node::TriplePreorderBifurcating(
    std::function<void(const Node*, const Node*, const Node*)> f) const {
  if (IsLeaf()) {
    return;
  }  // else
  std::stack<std::pair<const Node*, bool>> stack;
  stack.push({this, false});
  const Node* node;
  bool visited;
  while (stack.size()) {
    // Here we visit each node twice, once for each orientation.
    std::tie(node, visited) = stack.top();
    stack.pop();
    const auto& children = node->Children();
    Assert(children.size() == 2,
           "TriplePreorderBifurcating expects a bifurcating tree.");
    if (visited) {
      // We've already visited this node once, so do the second orientation.
      f(children[1].get(), children[0].get(), node);
      // Next traverse the right child.
      if (!children[1]->IsLeaf()) {
        stack.push({children[1].get(), false});
      }
    } else {
      // We are visiting this node for the first time.
      // Apply f in the first orientation.
      f(children[0].get(), children[1].get(), node);
      // Then set it up so it gets executed in the second orientation...
      stack.push({node, true});
      // ... after first traversing the left child.
      if (!children[0]->IsLeaf()) {
        stack.push({children[0].get(), false});
      }
    }
  }
}

// See the typedef of UnrootedPCSPFun to understand the argument type to this
// function, and `doc/svg/pcsp.svg` for a diagram that will greatly help you
// understand the implementation.
void Node::UnrootedPCSPPreorder(UnrootedPCSPFun f) const {
  this->TriplePreorder(
      // f_root
      [&f](const Node* node0, const Node* node1, const Node* node2) {
        // Virtual root on node2's edge, with subsplit pointing up.
        f(node2, false, node2, true, node0, false, node1, false, nullptr);
        if (!node2->IsLeaf()) {
          Assert(node2->Children().size() == 2,
                 "PCSPPreorder expects a bifurcating tree.");
          auto child0 = node2->Children()[0].get();
          auto child1 = node2->Children()[1].get();
          // Virtual root in node1.
          f(node0, false, node2, false, child0, false, child1, false, node1);
          // Virtual root in node0.
          f(node1, false, node2, false, child0, false, child1, false, node0);
          // Virtual root on node2's edge, with subsplit pointing down.
          f(node2, true, node2, false, child0, false, child1, false, nullptr);
          // Virtual root in child0.
          f(child1, false, node2, true, node0, false, node1, false, child0);
          // Virtual root in child1.
          f(child0, false, node2, true, node0, false, node1, false, child1);
        }
      },
      // f_internal
      [&f, this](const Node* node, const Node* sister, const Node* parent) {
        // Virtual root on node's edge, with subsplit pointing up.
        f(node, false, node, true, parent, true, sister, false, nullptr);
        if (!node->IsLeaf()) {
          Assert(node->Children().size() == 2,
                 "PCSPPreorder expects a bifurcating tree.");
          auto child0 = node->Children()[0].get();
          auto child1 = node->Children()[1].get();
          // Virtual root up the tree.
          f(sister, false, node, false, child0, false, child1, false, this);
          // Virtual root in sister.
          f(parent, true, node, false, child0, false, child1, false, sister);
          // Virtual root on node's edge, with subsplit pointing down.
          f(node, true, node, false, child0, false, child1, false, nullptr);
          // Virtual root in child0.
          f(child1, false, node, true, sister, false, parent, true, child0);
          // Virtual root in child1.
          f(child0, false, node, true, sister, false, parent, true, child1);
        }
      });
}

void Node::RootedPCSPPreorder(RootedPCSPFun f, bool allow_leaves) const {
  this->TriplePreorderBifurcating(
      [&f, &allow_leaves](const Node* node, const Node* sister, const Node* parent) {
        if (node->IsLeaf() && allow_leaves) {
          f(sister, node, nullptr, nullptr);
        } else if (!node->IsLeaf()) {
          Assert(node->Children().size() == 2,
                 "RootedPCSPPreorder expects a bifurcating tree.");
          auto child0 = node->Children()[0].get();
          auto child1 = node->Children()[1].get();
          f(sister, node, child0, child1);
        }
      });
}

void Node::RootedSisterAndLeafTraversal(TwoNodeFun f) const {
  this->TriplePreorderBifurcating(
      [&f](const Node* node, const Node* sister, const Node* parent) {
        if (node->IsLeaf()) {
          f(sister, node);
        }
      });
}

// This function assigns ids to the nodes of the topology: the leaves get
// their fixed ids (which we assume are contiguously numbered from 0 through
// the leaf count -1) and the rest get ordered according to a postorder
// traversal. Thus if the tree is bifurcating the root always has id equal to
// the number of nodes in the tree.
//
// This function returns a map that maps the tags to their ids.
TagSizeMap Node::Polish() {
  TagSizeMap tag_id_map;
  const size_t leaf_count = MaxLeafID() + 1;
  size_t next_id = leaf_count;
  MutablePostorder([&tag_id_map, &next_id, &leaf_count](Node* node) {
    if (node->IsLeaf()) {
      node->id_ = node->MaxLeafID();
      node->leaves_ = Bitset::Singleton(leaf_count, node->id_);
    } else {
      node->id_ = next_id;
      next_id++;
      node->leaves_ = Node::LeavesOf(node->Children());
    }
    SafeInsert(tag_id_map, node->Tag(), node->id_);
  });
  return tag_id_map;
}

std::string Node::Newick(std::function<std::string(const Node*)> node_labeler,
                         const DoubleVectorOption& branch_lengths) const {
  return NewickAux(node_labeler, branch_lengths) + ";";
}

std::string Node::NewickAux(std::function<std::string(const Node*)> node_labeler,
                            const DoubleVectorOption& branch_lengths) const {
  std::string str;
  if (IsLeaf()) {
    str.assign(node_labeler(this));
  } else {
    str.assign("(");
    for (auto iter = children_.begin(); iter != children_.end(); iter++) {
      if (iter != children_.begin()) {
        str.append(",");
      }
      str.append((*iter)->NewickAux(node_labeler, branch_lengths));
    }
    str.append(")");
    str.append(node_labeler(this));
  }
  if (branch_lengths) {
    Assert(Id() < (*branch_lengths).size(),
           "branch_lengths vector is of insufficient length in NewickAux.");
    // ostringstream is the way to get scientific notation using the STL.
    std::ostringstream str_stream;
    str_stream << (*branch_lengths)[Id()];
    str.append(":" + str_stream.str());
  }
  return str;
}

std::string Node::Newick(const DoubleVectorOption& branch_lengths,
                         const TagStringMapOption& node_labels, bool show_tags) const {
  return Newick(
      [&node_labels, &show_tags](const Node* node) {
        if (node->IsLeaf()) {
          if (node_labels) {
            return (*node_labels).at(node->Tag());
          } else if (show_tags) {
            return node->TagString();
          } else {
            return std::to_string(node->MaxLeafID());
          }
        } else {
          if (show_tags) {
            return node->TagString();
          }
        }
        return std::string("");
      },
      branch_lengths);
}

std::vector<size_t> Node::ParentIdVector() const {
  std::vector<size_t> ids(Id());
  Postorder([&ids](const Node* node) {
    if (!node->IsLeaf()) {
      for (const auto& child : node->Children()) {
        Assert(child->Id() < ids.size(), "Problematic ids in ParentIdVector.");
        ids[child->Id()] = node->Id();
      }
    }
  });
  return ids;
}

Node::NodePtr Node::Deroot() {
  Assert(LeafCount() >= 3, "Node::Deroot expects a tree with at least 3 tips.");
  Assert(Children().size() == 2, "Can't deroot a non-bifurcating tree.");
  auto deroot = [](const NodePtr other_child, const NodePtr has_descendants) {
    // Make a vector copy by passing a vector in.
    NodePtrVec children(has_descendants->Children());
    children.push_back(other_child);
    // has_descendants' id is now available.
    return Join(children, has_descendants->Id());
  };
  if (children_[1]->LeafCount() == 1) {
    return deroot(children_[1], children_[0]);
  }  // else
  return deroot(children_[0], children_[1]);
}

Bitset Node::LeavesOf(const Node::NodePtrVec& children) {
  Assert(children.size() > 0, "Need children in Node::LeavesOf.");
  Bitset leaves(children[0]->Leaves());
  for (size_t i = 1; i < children.size(); i++) {
    leaves |= children[i]->Leaves();
  }
  return leaves;
}

// Class methods
Node::NodePtr Node::Leaf(uint32_t id, Bitset leaves) {
  return std::make_shared<Node>(id, leaves);
}
Node::NodePtr Node::Leaf(uint32_t id, size_t taxon_count) {
  return Node::Leaf(id, Bitset::Singleton(taxon_count, id));
}
Node::NodePtr Node::Join(NodePtrVec children, size_t id) {
  return std::make_shared<Node>(children, id, Node::LeavesOf(children));
}
Node::NodePtr Node::Join(NodePtr left, NodePtr right, size_t id) {
  return Join(std::vector<NodePtr>({left, right}), id);
}

Node::NodePtr Node::OfParentIdVector(const std::vector<size_t>& ids) {
  // We will fill this map with the ids of the descendants.
  std::unordered_map<size_t, std::vector<size_t>> downward_ids;
  for (size_t child_id = 0; child_id < ids.size(); child_id++) {
    const auto& parent_id = ids[child_id];
    auto search = downward_ids.find(parent_id);
    if (search == downward_ids.end()) {
      // The first time we have seen this parent.
      std::vector<size_t> child_ids({child_id});
      SafeInsert(downward_ids, parent_id, std::move(child_ids));
    } else {
      // We've seen the parent before, so append the child to the parent's
      // vector of descendants.
      search->second.push_back(child_id);
    }
  }
  // The leaf count is equal to the smallest non-leaf index, i.e. a parent
  // index.
  size_t leaf_count = *std::min_element(ids.begin(), ids.end());
  std::function<NodePtr(size_t)> build_tree = [&build_tree, &downward_ids,
                                               leaf_count](size_t current_id) {
    auto search = downward_ids.find(current_id);
    if (search == downward_ids.end()) {
      // We assume that anything not in the map is a leaf, because leaves
      // don't have any children.
      return Leaf(static_cast<uint32_t>(current_id),
                  Bitset::Singleton(leaf_count, current_id));
    } else {
      const auto& children_ids = search->second;
      std::vector<NodePtr> children;
      for (const auto& child_id : children_ids) {
        children.push_back(build_tree(child_id));
      }
      return Join(children, current_id);
    }
  };
  // We assume that the maximum id of the tree is the length of the input
  // id array. That makes sense because the root does not have a parent, so
  // is the first "missing" entry in the input id array.
  return build_tree(ids.size());
}

Node::NodePtrVec Node::ExampleTopologies() {
  NodePtrVec topologies = {
      // 0: (0,1,(2,3))
      Join(std::vector<NodePtr>({Leaf(0), Leaf(1), Join(Leaf(2), Leaf(3))})),
      // 1; (0,1,(2,3)) again
      Join(std::vector<NodePtr>({Leaf(1), Leaf(0), Join(Leaf(3), Leaf(2))})),
      // 2: (0,2,(1,3))
      Join(std::vector<NodePtr>({Leaf(0), Leaf(2), Join(Leaf(1), Leaf(3))})),
      // 3: (0,(1,(2,3)))
      Join(std::vector<NodePtr>({Leaf(0), Join(Leaf(1), Join(Leaf(2), Leaf(3)))})),
      // 4: ((0,(2,3)),1)
      Join(std::vector<NodePtr>({Join(Leaf(0), Join(Leaf(2), Leaf(3))), Leaf(1)}))};
  for (auto& topology : topologies) {
    topology->Polish();
  }
  return topologies;
}

Node::NodePtr Node::Ladder(uint32_t leaf_count) {
  Assert(leaf_count > 0, "leaf_count should be positive in Node::Ladder.");
  NodePtr node = Leaf(0);
  for (uint32_t i = 1; i < leaf_count; i++) {
    node = Join(Leaf(i), node);
  }
  node->Polish();
  return node;
}

inline uint32_t Node::SOHash(uint32_t x) {
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = (x >> 16) ^ x;
  return x;
}

inline size_t Node::SORotate(size_t n, uint32_t c) {
  const uint32_t mask = (CHAR_BIT * sizeof(n) - 1);  // assumes width is a power of 2.
  // assert ( (c<=mask) &&"rotate by type width or more");
  c &= mask;
  return (n << c) | (n >> ((-c) & mask));
}

SizeVectorVector Node::IdsAbove() const {
  SizeVectorVector ids_above(Id() + 1);
  SizeVector mutable_ids;
  DepthFirst(
      [&ids_above, &mutable_ids](const Node* node) {
        // Store the current set of ids above.
        ids_above[node->Id()] = SizeVector(mutable_ids);
        // As we travel down the tree, the current node will be above.
        mutable_ids.push_back(node->Id());
      },
      // Going back up the tree, so remove the current node's id.
      [&mutable_ids](const Node*) { mutable_ids.pop_back(); });
  return ids_above;
}

std::string Node::NodeIdAndLeavesToString() const {
  std::stringstream os;
  os << "{ id: " << Id();
  os << ", leaves: " << Leaves();
  SizeVector child_ids;
  for (const auto& child : children_) {
    child_ids.push_back(child.get()->Id());
  }
  os << ", children: " << child_ids << " }";
  return os.str();
}

std::string Node::NodeIdAndLeavesToStringForTopology() const {
  std::stringstream os;
  os << "[ " << std::endl;

  Preorder(
      // pre function
      [&os](const Node* node_ptr) {
        os << "\t" << node_ptr->NodeIdAndLeavesToString() << ", " << std::endl;
      });

  os << " ]" << std::endl;
  return os.str();
}
