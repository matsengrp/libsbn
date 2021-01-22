// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// A "tidy" subsplit DAG has a notion of clean and dirty vectors.

#ifndef SRC_TIDY_SUBSPLIT_DAG_HPP_
#define SRC_TIDY_SUBSPLIT_DAG_HPP_

#include "subsplit_dag.hpp"

class TidySubsplitDAG : public SubsplitDAG {
 public:
  TidySubsplitDAG();
  // This constructor is really just meant for testing.
  explicit TidySubsplitDAG(EigenMatrixXb above);
  explicit TidySubsplitDAG(const RootedTreeCollection &tree_collection);

  EigenArrayXbRef BelowNode(size_t node_idx);
  EigenArrayXb AboveNode(size_t node_idx);
  void JoinBelow(size_t dst, size_t src1, size_t src2);

  void PrintBelowMatrix();

 private:
  // above_.(i,j) is true if i is above j, otherwise it's false.
  EigenMatrixXb above_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("TidySubsplitDAG: slicing") {
  EigenMatrixXb above(5, 5);
  above.setIdentity();
  auto dag = TidySubsplitDAG(above);

  dag.JoinBelow(1, 0, 2);
  dag.JoinBelow(3, 1, 4);

  dag.PrintBelowMatrix();
  // CHECK_EQ(GenericToString(dag.AboveNode(1)), "[0, 1, 0]\n");
  // CHECK_EQ(GenericToString(dag.BelowNode(1)), "[0, 1, 1]\n");
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_TIDY_SUBSPLIT_DAG_HPP_