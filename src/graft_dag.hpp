// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// The GraftDAG is a proposed addition (graft) to SubsplitDAG (host), which we
// can perform computations on without the need for adding nodes and edges and
// reindexing the full DAG.

#include "gp_dag.hpp"
#include "subsplit_dag.hpp"

#pragma once

class GraftDAG : public SubsplitDAG {
 public:
  // ** Constructors:

  // Initialize empty GraftDAG.
  GraftDAG(SubsplitDAG &dag);

  // ** Comparators

  // Uses same method of comparison as SubsplitDAG (node and edge sets).
  int Compare(const GraftDAG &other) const;
  static int Compare(const GraftDAG &lhs, const GraftDAG &rhs);
  // Treats GraftDAG as completed DAG to compare against normal SubsplitDAG.
  int CompareToDAG(const SubsplitDAG &other) const;
  static int CompareToDAG(const GraftDAG &lhs, const SubsplitDAG &rhs);

  // ** Modify GraftDAG

  // Clear all nodes and edges from graft for reuse.
  void RemoveAllGrafts();

  // ** Getters

  // Get pointer to the host DAG.
  const SubsplitDAG &GetHostDAG() const;

  // ** Counts

  // Total number of nodes in graft only.
  size_t GraftNodeCount() const;
  // Total number of nodes in host DAG only.
  size_t HostNodeCount() const;
  // Total number of edges in graft only.
  size_t GraftEdgeCount() const;
  // Total number of edges in host DAG only.
  size_t HostEdgeCount() const;
  // Check if a node is from the host, otherwise from the graft.
  bool IsNodeFromHost(size_t node_id) const;
  // Check if an edge is from the host, otherwise from the graft.
  bool IsEdgeFromHost(size_t edge_id) const;

  // ** Contains

  // Checks whether the node is in the graft only.
  bool ContainsGraftNode(const Bitset node_subsplit) const;
  bool ContainsGraftNode(const size_t node_id) const;
  // Checks whether the edge is in the graft only.
  bool ContainsGraftEdge(const size_t parent_id, const size_t child_id) const;

 protected:
  // DAG that the graft is proposed to be connected to.
  SubsplitDAG &host_dag_;
};
