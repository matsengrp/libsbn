// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_UNROOTED_TREE_COLLECTION_HPP_
#define SRC_UNROOTED_TREE_COLLECTION_HPP_

#include "tree_collection.hpp"
#include "unrooted_tree.hpp"

using PreUnrootedTreeCollection = GenericTreeCollection<UnrootedTree>;

class UnrootedTreeCollection : public PreUnrootedTreeCollection {
 public:
  // Inherit all constructors.
  using PreUnrootedTreeCollection::PreUnrootedTreeCollection;

  static UnrootedTreeCollection OfTreeCollection(const TreeCollection& trees);
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("UnrootedTreeCollection") {}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_UNROOTED_TREE_COLLECTION_HPP_
