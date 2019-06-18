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


class Node {

 typedef std::shared_ptr<Node> NodePtr;
 typedef std::vector<NodePtr> NodePtrVec;

 private:
  NodePtrVec children_;
  int id_;


 public:
  Node(int id) : children_({}), id_(id) {}

  Node(NodePtrVec children, int id) : children_(children), id_(id) {}

  Node(Node & left, Node & right, int id) {
    auto leftptr = std::make_shared<Node>(left);
    auto rightptr = std::make_shared<Node>(right);
    children_ = {leftptr, rightptr};
    id_ = id;
  }

  Node(NodePtr left, NodePtr right, int id) {
    children_ = {left, right};
    id_ = id;
  }

  ~Node() { std::cout << "Destroying node " << id_ << std::endl; }

  int get_id() { return id_; }
  bool is_leaf() { return children_.empty(); }

  int n_leaves() {
    if (is_leaf()) {
      return 1;
    }

    int count = 0;
    for (auto child : children_) {
      count += child->n_leaves();
    }

    return count;
  }
};


TEST_CASE("Trying out Node") {
  auto l1 = std::make_shared<Node>(1);
  auto l2 = std::make_shared<Node>(2);
  auto t = std::make_shared<Node>(l1, l2, 3);

  CHECK(t->n_leaves() == 2);
}

#endif
