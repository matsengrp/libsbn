// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_ALIGNMENT_HPP_
#define SRC_ALIGNMENT_HPP_

#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include "typedefs.hpp"

// The mutable nature of this class is a little ugly.
// However, its only purpose is to sit in a libsbn Instance, which is all about
// mutable state.
class Alignment {
 public:
  Alignment() {}
  explicit Alignment(StringStringMap data) : data_(data) {}

  StringStringMap Data() const { return data_; }
  size_t SequenceCount() const { return data_.size(); }
  size_t Length() const {
    assert(SequenceCount() > 0);
    return data_.begin()->second.size();
  }

  // Is the alignment non-empty and do all sequences have the same length?
  bool IsValid() const {
    if (data_.size() == 0) {
      return false;
    }
    size_t length = data_.begin()->second.size();
    for (auto iter = data_.begin(); iter != data_.end(); ++iter) {
      if (length != iter->second.size()) {
        return false;
      }
    }
    return true;
  }

  std::string at(const std::string &taxon) const {
    auto search = data_.find(taxon);
    if (search != data_.end()) {
      return search->second;
    } else {
      std::cerr << "Taxon '" << search->first << "' not found in alignment.\n";
      abort();
    }
  }

  // An edited version of
  // https://stackoverflow.com/questions/35251635/fasta-reader-written-in-c
  // which seems like it was originally taken from
  // http://rosettacode.org/wiki/FASTA_format#C.2B.2B
  void ReadFasta(std::string fname) {
    StringStringMap &data = this->data_;
    data.clear();
    auto insert = [&data](std::string taxon, std::string sequence) {
      if (!taxon.empty()) {
        assert(data.insert({taxon, sequence}).second);
      }
    };
    std::ifstream input(fname);
    if (!input.good()) {
      std::cerr << "Could not open '" << fname << "'\n";
      abort();
    }
    std::string line, taxon, sequence;
    while (std::getline(input, line)) {
      if (line.empty()) continue;
      if (line[0] == '>') {
        insert(taxon, sequence);
        taxon = line.substr(1);
        sequence.clear();
      } else {
        sequence += line;
      }
    }
    // Insert the last taxon, sequence pair.
    insert(taxon, sequence);
    if (!IsValid()) {
      std::cerr << "Sequences of the alignment are not all the same length.\n";
      abort();
    }
  }

 private:
  StringStringMap data_;
};

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("Alignment") {
  Alignment alignment;
  alignment.ReadFasta("data/hello.fasta");
  StringStringMap correct({{"mars", "CCGAG-AGCAGCAATGGAT-GAGGCATGGCG"},
                           {"saturn", "GCGCGCAGCTGCTGTAGATGGAGGCATGACG"},
                           {"jupiter", "GCGCGCAGCAGCTGTGGATGGAAGGATGACG"}});
  CHECK_EQ(correct, alignment.Data());
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif  // SRC_ALIGNMENT_HPP_

