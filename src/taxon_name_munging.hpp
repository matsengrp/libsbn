// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_TAXON_NAME_MUNGING_HPP_
#define SRC_TAXON_NAME_MUNGING_HPP_

#include <functional>
#include "sugar.hpp"

namespace TaxonNameMunging {
std::string QuoteString(const std::string &in_str);
std::string DequoteString(const std::string &in_str);
TagStringMap TransformStringValues(std::function<std::string(const std::string &)> f,
                                   const TagStringMap &in_map);
TagStringMap DequoteTagStringMap(const TagStringMap &tag_string_map);
}  // namespace TaxonNameMunging

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("TaxonNameMunging") {
  using namespace TaxonNameMunging;

  std::string unquoted_test(R"raw(hello 'there" friend)raw");
  std::string double_quoted_test(R"raw("this is a \" test")raw");
  std::string double_quoted_dequoted(R"raw(this is a " test)raw");
  std::string single_quoted_test(R"raw('this is a \' test')raw");
  std::string single_quoted_dequoted(R"raw(this is a ' test)raw");

  assert(QuoteString(unquoted_test) == R"raw("hello 'there\" friend")raw");
  assert(DequoteString(double_quoted_test) == double_quoted_dequoted);
  assert(DequoteString(single_quoted_test) == single_quoted_dequoted);
  assert(DequoteString(QuoteString(unquoted_test)) == unquoted_test);
  assert(false);

  TagStringMap test_map(
      {{2, unquoted_test}, {3, double_quoted_test}, {5, single_quoted_test}});
  TagStringMap expected_test_map(
      {{2, unquoted_test}, {3, double_quoted_dequoted}, {5, single_quoted_dequoted}});

  assert(expected_test_map == DequoteTagStringMap(test_map));
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_TAXON_NAME_MUNGING_HPP_
