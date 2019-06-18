#ifndef __SBN_HPP
#define __SBN_HPP

#include "doctest.h"

#include <iostream>
#include <memory>
#include <vector>

class MyClass {
 private:
  int m_i;

 public:
  void int_set(int i);

  int int_get();
};


// Class for a tree.
// Has nodes with unsigned integer ids.
// These ids have to increase as we go towards the root.
class Node {
  typedef std::shared_ptr<Node> NodePtr;
  typedef std::vector<NodePtr> NodePtrVec;

 private:
  NodePtrVec children_;
  // TODO: "NodeId" type rather than unsigned int?
  unsigned int id_;

  // Make copy constructors private to eliminate copying.
  Node(const Node&);
  Node& operator=(const Node&);


 public:
  Node(unsigned int id) : children_({}), id_(id) {}
  Node(NodePtrVec children, unsigned int id) : children_(children), id_(id) {
    REQUIRE_MESSAGE(MaxChildIdx(children) < id,
                    "Nodes must have a larger index than their children.");
  }
  Node(NodePtr left, NodePtr right, unsigned int id)
      : Node({left, right}, id) {}
  ~Node() { std::cout << "Destroying node " << id_ << std::endl; }

  unsigned int GetId() { return id_; }
  bool IsLeaf() { return children_.empty(); }

  void PreOrder(std::function<void(Node*)> f) {
    f(this);
    for (auto child : children_) {
      child->PreOrder(f);
    }
  }

  void PostOrder(std::function<void(Node*)> f) {
    for (auto child : children_) {
      child->PostOrder(f);
    }
    f(this);
  }

  unsigned int LeafCount() {
    unsigned int count = 0;
    PreOrder([&count](Node* node) { count += node->IsLeaf(); });
    return count;
  }

  // Class methods
  static NodePtr Leaf(int id) { return std::make_shared<Node>(id); }
  static NodePtr Join(NodePtr left, NodePtr right, int id) {
    return std::make_shared<Node>(left, right, id);
  };
  static unsigned int MaxChildIdx(NodePtrVec children) {
    REQUIRE(~children.empty());
    // 0 is the smallest value for an unsigned integer.
    unsigned int result = 0;
    for (auto child : children) {
      result = std::max(child->GetId(), result);
    }
    return result;
  }
};


TEST_CASE("Trying out Node") {
  auto l0 = Node::Leaf(0);
  auto l1 = Node::Leaf(1);
  auto t = Node::Join(l0, l1, 2);

  // Should fail.
  // Node::Join(l1, l2, 0);

  REQUIRE(t->LeafCount() == 2);

  auto print_pos = [](Node* t) {
    std::cout << "I'm at " << t->GetId() << std::endl;
  };
  t->PreOrder(print_pos);
  t->PostOrder(print_pos);

}

#endif
