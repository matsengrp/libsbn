// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_LIBSBN_HPP_
#define SRC_LIBSBN_HPP_

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <cmath>
#include <queue>
#include <random>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>
#include "alignment.hpp"
#include "beagle.hpp"
#include "build.hpp"
#include "driver.hpp"
#include "task_processor.hpp"
#include "tree.hpp"

namespace py = pybind11;

typedef std::unordered_map<std::string, float> StringFloatMap;
typedef std::unordered_map<std::string, uint32_t> StringUInt32Map;
typedef std::unordered_map<std::string, std::pair<uint32_t, uint32_t>>
    StringUInt32PairMap;

template <typename T>
StringUInt32Map StringUInt32MapOf(T m) {
  StringUInt32Map m_str;
  for (const auto &iter : m) {
    m_str[iter.first.ToString()] = iter.second;
  }
  return m_str;
}

struct SBNInstance {
  std::string name_;
  TreeCollection::TreeCollectionPtr tree_collection_;
  Alignment alignment_;
  CharIntMap symbol_table_;
  std::vector<beagle::BeagleInstance> beagle_instances_;
  std::vector<double> sbn_probs_;
  BitsetUInt32Map indexer_;
  UInt32BitsetMap reverse_indexer_;
  UInt32BitsetMap index_to_child_;
  BitsetUInt32PairMap range_indexer_;
  // The first index after the rootsplit block.
  size_t rootsplit_index_end_;
  static std::random_device random_device_;
  static std::mt19937 random_generator_;

  // ** Initialization, destruction, and status
  explicit SBNInstance(const std::string &name)
      : name_(name), symbol_table_(beagle::GetSymbolTable()) {}

  ~SBNInstance() { FinalizeBeagleInstances(); }

  // Finalize means to release memory.
  void FinalizeBeagleInstances() {
    for (const auto &beagle_instance : beagle_instances_) {
      assert(beagleFinalizeInstance(beagle_instance) == 0);
    }
    beagle_instances_.clear();
  }

  size_t TreeCount() const { return tree_collection_->TreeCount(); }
  void PrintStatus() {
    std::cout << "Status for instance '" << name_ << "':\n";
    if (tree_collection_) {
      std::cout << TreeCount() << " unique tree topologies loaded on "
                << tree_collection_->TaxonCount() << " leaves.\n";
    } else {
      std::cout << "No trees loaded.\n";
    }
    std::cout << alignment_.Data().size() << " sequences loaded.\n";
  }

  // ** Building SBN-related items

  void ProcessLoadedTrees() {
    uint32_t index = 0;
    auto counter = tree_collection_->TopologyCounter();
    indexer_.clear();
    reverse_indexer_.clear();
    range_indexer_.clear();
    index_to_child_.clear();
    // Start by adding the rootsplits.
    for (const auto &iter : RootsplitCounterOf(counter)) {
      indexer_[iter.first] = index;
      index++;
    }
    rootsplit_index_end_ = index;
    // Now we make a map that maps parent subsplit to the collection of PCSSs
    // that contain it.
    // TODO think about move and make consistent below
    BitsetBitsetVectorMap parent_to_pcss_vector_map;
    for (const auto &iter : PCSSCounterOf(counter)) {
      const auto &pcss = iter.first;
      auto parent = pcss.PCSSParent();
      auto search = parent_to_pcss_vector_map.find(parent);
      if (search == parent_to_pcss_vector_map.end()) {
        std::vector<Bitset> pcss_singleton({pcss});
        assert(parent_to_pcss_vector_map
                   .insert({std::move(parent), std::move(pcss_singleton)})
                   .second);
      } else {
        search->second.push_back(std::move(pcss));
      }
    }
    for (auto iter : parent_to_pcss_vector_map) {
      auto parent = iter.first;
      auto pcss_vector = iter.second;
      // The range indexer maps the parent to the range of values used by that
      // parent.
      assert(
          range_indexer_.insert({parent, {index, index + pcss_vector.size()}})
              .second);
      // We insert the corresponding PCSSs into the indexer.
      for (const auto &pcss : pcss_vector) {
        assert(indexer_.insert({pcss, index}).second);
        index++;
      }
    }
    sbn_probs_ = std::vector<double>(indexer_.size(), 1.);
    for (const auto &iter : indexer_) {
      assert(reverse_indexer_.insert({iter.second, iter.first}).second);
      assert(reverse_indexer_.insert({iter.second, iter.first.PCSSChild()})
                 .second);
    }
  }

  uint32_t SampleIndex(std::pair<uint32_t, uint32_t> range) const {
    assert(range.first < range.second);
    assert(range.second <= sbn_probs_.size());
    std::discrete_distribution<> distribution(
        sbn_probs_.begin() + range.first, sbn_probs_.begin() + range.second);
    return static_cast<uint32_t>(distribution(random_generator_));
  }

  // TODO rootsplit naming.
  // This function samples a tree by first sampling the rootsplit, and then
  // calling the recursive form of SampleTree.
  Node::NodePtr SampleTree() const {
    auto rootsplit_index =
        SampleIndex(std::pair<uint32_t, uint32_t>(0, rootsplit_index_end_));
    const Bitset &rootsplit = reverse_indexer_.at(rootsplit_index);
    return SampleTree(rootsplit);
  }

  Node::NodePtr SampleTree(const Bitset &parent_split) const {
    auto process_split = [this](const Bitset &split) {
      auto singleton_option = split.SplitChunk(1).SingletonOption();
      if (singleton_option) {
        return Node::Leaf(*singleton_option);
      }  // else
      auto child_split =
          index_to_child_.at(SampleIndex(range_indexer_.at(split)));
      return SampleTree(child_split);
    };
    std::cout << parent_split.ToString() << std::endl;
    return Node::Join(process_split(parent_split),
                      process_split(parent_split.SisterExchange()));
  }

  // TODO(erick) replace with something interesting.
  double SBNTotalProb() {
    double total = 0;
    for (const auto &prob : sbn_probs_) {
      total += prob;
    }
    return total;
  }

  // ** I/O

  std::tuple<StringUInt32Map, StringUInt32PairMap> GetIndexers() {
    auto indexer_str = StringUInt32MapOf(indexer_);
    StringUInt32PairMap range_indexer_str;
    for (const auto &iter : range_indexer_) {
      assert(range_indexer_str.insert({iter.first.ToString(), iter.second})
                 .second);
    }
    assert(range_indexer_str.insert({"rootsplit", {0, rootsplit_index_end_}})
               .second);
    return std::tie(indexer_str, range_indexer_str);
  }

  // This function is really just for testing-- it recomputes from scratch.
  std::pair<StringUInt32Map, StringUInt32Map> SplitCounters() {
    auto counter = tree_collection_->TopologyCounter();
    return {StringUInt32MapOf(RootsplitCounterOf(counter)),
            StringUInt32MapOf(PCSSCounterOf(counter))};
  }

  void ReadNewickFile(std::string fname) {
    Driver driver;
    tree_collection_ = driver.ParseNewickFile(fname);
  }

  void ReadNexusFile(std::string fname) {
    Driver driver;
    tree_collection_ = driver.ParseNexusFile(fname);
  }

  void ReadFastaFile(std::string fname) { alignment_.ReadFasta(fname); }

  // ** Phylogenetic likelihood

  void MakeBeagleInstances(int instance_count) {
    // Start by clearing out any existing instances.
    FinalizeBeagleInstances();
    if (alignment_.SequenceCount() == 0) {
      std::cerr << "Load an alignment into your SBNInstance on which you wish "
                   "to calculate phylogenetic likelihoods.\n";
      abort();
    }
    if (TreeCount() == 0) {
      std::cerr << "Load some trees into your SBNInstance on which you wish to "
                   "calculate phylogenetic likelihoods.\n";
      abort();
    }
    for (auto i = 0; i < instance_count; i++) {
      auto beagle_instance = beagle::CreateInstance(alignment_);
      beagle::SetJCModel(beagle_instance);
      beagle_instances_.push_back(beagle_instance);
      beagle::PrepareBeagleInstance(beagle_instance, tree_collection_,
                                    alignment_, symbol_table_);
    }
  }

  std::vector<double> TreeLogLikelihoods() {
    if (beagle_instances_.size() == 0) {
      std::cerr << "Please add some BEAGLE instances that can be used for "
                   "computation.\n";
      abort();
    }
    std::vector<double> results(tree_collection_->TreeCount());
    std::queue<beagle::BeagleInstance> instance_queue;
    for (auto instance : beagle_instances_) {
      instance_queue.push(instance);
    }
    std::queue<size_t> tree_number_queue;
    for (size_t i = 0; i < tree_collection_->TreeCount(); i++) {
      tree_number_queue.push(i);
    }
    TaskProcessor<beagle::BeagleInstance, size_t> task_processor(
        instance_queue, tree_number_queue,
        [&results, &tree_collection = tree_collection_ ](
            beagle::BeagleInstance beagle_instance, size_t tree_number) {
          results[tree_number] = beagle::TreeLogLikelihood(
              tree_collection->GetTree(tree_number), beagle_instance);
        });
    return results;
  }
};

// Here we initialize our static random number generator.
std::random_device SBNInstance::random_device_;
std::mt19937 SBNInstance::random_generator_(SBNInstance::random_device_());

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("libsbn") {
  SBNInstance inst("charlie");
  inst.ReadNewickFile("data/five_taxon.nwk");
  // Reading one file after another checks that we've cleared out state.
  inst.ReadNewickFile("data/hello.nwk");
  inst.ReadFastaFile("data/hello.fasta");
  inst.MakeBeagleInstances(2);
  for (auto ll : inst.TreeLogLikelihoods()) {
    CHECK_LT(abs(ll - -84.852358), 0.000001);
  }
  // inst.ProcessLoadedTrees();
  // auto tree = inst.SampleTree();

  inst.ReadNexusFile("data/DS1.subsampled_10.t");
  inst.ReadFastaFile("data/DS1.fasta");
  inst.MakeBeagleInstances(2);
  auto likelihoods = inst.TreeLogLikelihoods();
  std::vector<double> pybeagle_likelihoods(
      {-14582.995273982739, -6911.294207416366, -6916.880235529542,
       -6904.016888831189, -6915.055570693576, -6915.50496696512,
       -6910.958836661867, -6909.02639968063, -6912.967861935749,
       -6910.7871105783515});
  for (size_t i = 0; i < likelihoods.size(); i++) {
    CHECK_LT(abs(likelihoods[i] - pybeagle_likelihoods[i]), 0.00011);
  }
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_LIBSBN_HPP_
