// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_LIBSBN_HPP_
#define SRC_LIBSBN_HPP_

#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <algorithm>
#include <cmath>
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
#include "prettyprint.hpp"
#include "tree.hpp"

namespace py = pybind11;

typedef std::unordered_map<std::string, float> StringFloatMap;
typedef std::unordered_map<std::string, uint32_t> StringUInt32Map;
typedef std::unordered_map<std::string, std::pair<uint32_t, uint32_t>>
    StringUInt32PairMap;
typedef std::unordered_map<uint32_t, std::string> UInt32StringMap;
typedef std::unordered_map<std::string,
                           std::unordered_map<std::string, uint32_t>>
    StringPCSSMap;
typedef std::vector<uint32_t> PCSSIndexVector;

template <typename T>
StringUInt32Map StringUInt32MapOf(T m) {
  StringUInt32Map m_str;
  for (const auto &iter : m) {
    m_str[iter.first.ToString()] = iter.second;
  }
  return m_str;
}

StringPCSSMap StringPCSSMapOf(PCSSDict d) {
  StringPCSSMap d_str;
  for (const auto &iter : d) {
    d_str[iter.first.ToString()] = StringUInt32MapOf(iter.second);
  }
  return d_str;
}

struct SBNInstance {
  std::string name_;
  // Things that get loaded in or sampled from SBNs.
  TreeCollection tree_collection_;
  Alignment alignment_;
  // Beagly bits.
  CharIntMap symbol_table_;
  std::vector<beagle::BeagleInstance> beagle_instances_;
  size_t beagle_leaf_count_;
  size_t beagle_site_count_;
  // A vector that contains all of the SBN-related probabilities.
  std::vector<double> sbn_probs_;
  // A map that indexes these probabilities: rootsplits are at the beginning,
  // and PCSS bitsets are at the end.
  // The collection of rootsplits, with the same indexing as in the indexer_.
  BitsetVector rootsplits_;
  // The first index after the rootsplit block in sbn_probs_.
  size_t rootsplit_index_end_;
  BitsetUInt32Map indexer_;
  // A map going from the index of a PCSS to its child.
  UInt32BitsetMap index_to_child_;
  // A map going from a parent subsplit to the range of indices in sbn_probs_
  // with its children.
  BitsetUInt32PairMap parent_to_range_;
  // Random bits.
  static std::random_device random_device_;
  static std::mt19937 random_generator_;

  // ** Initialization, destruction, and status
  explicit SBNInstance(const std::string &name)
      : name_(name),
        symbol_table_(beagle::GetSymbolTable()),
        beagle_leaf_count_(0),
        beagle_site_count_(0) {}

  ~SBNInstance() { FinalizeBeagleInstances(); }

  // Finalize means to release memory.
  void FinalizeBeagleInstances() {
    for (const auto &beagle_instance : beagle_instances_) {
      assert(beagleFinalizeInstance(beagle_instance) == 0);
    }
    beagle_instances_.clear();
    beagle_leaf_count_ = 0;
    beagle_site_count_ = 0;
  }

  size_t TreeCount() const { return tree_collection_.TreeCount(); }
  void PrintStatus() {
    std::cout << "Status for instance '" << name_ << "':\n";
    if (tree_collection_.TreeCount()) {
      std::cout << TreeCount() << " unique tree topologies loaded on "
                << tree_collection_.TaxonCount() << " leaves.\n";
    } else {
      std::cout << "No trees loaded.\n";
    }
    std::cout << alignment_.Data().size() << " sequences loaded.\n";
  }

  // ** Building SBN-related items

  void ProcessLoadedTrees() {
    uint32_t index = 0;
    auto counter = tree_collection_.TopologyCounter();
    // See above for the definitions of these members.
    sbn_probs_.clear();
    rootsplits_.clear();
    indexer_.clear();
    index_to_child_.clear();
    parent_to_range_.clear();
    // Start by adding the rootsplits.
    for (const auto &iter : RootsplitCounterOf(counter)) {
      assert(indexer_.insert({iter.first, index}).second);
      rootsplits_.push_back(iter.first);
      index++;
    }
    rootsplit_index_end_ = index;
    // Now add the PCSSs.
    for (const auto &iter : PCSSCounterOf(counter)) {
      const auto &parent = iter.first;
      const auto &child_counter = iter.second;
      assert(parent_to_range_
                 .insert({parent, {index, index + child_counter.size()}})
                 .second);
      for (const auto &child_iter : child_counter) {
        const auto &child = child_iter.first;
        assert(indexer_.insert({parent + child, index}).second);
        assert(index_to_child_
                   .insert({index, Bitset::ChildSubsplit(parent, child)})
                   .second);
        index++;
      }
    }
    sbn_probs_ = std::vector<double>(index, 1.);
  }

  void PrintSupports() {
    std::vector<std::string> to_print(indexer_.size());
    for (const auto &iter : indexer_) {
      if (iter.second < rootsplit_index_end_) {
        to_print[iter.second] = iter.first.ToString();
      } else {
        to_print[iter.second] = iter.first.PCSSToString();
      }
    }
    for (size_t i = 0; i < to_print.size(); i++) {
      std::cout << i << "\t" << to_print[i] << std::endl;
    }
  }

  // Sample an integer index in [range.first, range.second) according to
  // sbn_probs_.
  uint32_t SampleIndex(std::pair<uint32_t, uint32_t> range) const {
    assert(range.first < range.second);
    assert(range.second <= sbn_probs_.size());
    std::discrete_distribution<> distribution(
        sbn_probs_.begin() + range.first, sbn_probs_.begin() + range.second);
    // We have to add on range.first because we have taken a slice of the full
    // array, and the sampler treats the beginning of this slice as zero.
    auto result =
        range.first + static_cast<uint32_t>(distribution(random_generator_));
    assert(result < range.second);
    return result;
  }

  // This function samples a tree by first sampling the rootsplit, and then
  // calling the recursive form of SampleTopology.
  Node::NodePtr SampleTopology() const {
    // Start by sampling a rootsplit.
    uint32_t rootsplit_index =
        SampleIndex(std::pair<uint32_t, uint32_t>(0, rootsplit_index_end_));
    const Bitset &rootsplit = rootsplits_.at(rootsplit_index);
    // The addition below turns the rootsplit into a subsplit.
    auto topology = SampleTopology(rootsplit + ~rootsplit)->Deroot();
    topology->Reindex();
    return topology;
  }

  // The input to this function is a parent subsplit (of length 2n).
  Node::NodePtr SampleTopology(const Bitset &parent_subsplit) const {
    auto process_subsplit = [this](const Bitset &parent) {
      auto singleton_option = parent.SplitChunk(1).SingletonOption();
      if (singleton_option) {
        return Node::Leaf(*singleton_option);
      }  // else
      auto child_index = SampleIndex(parent_to_range_.at(parent));
      return SampleTopology(index_to_child_.at(child_index));
    };
    return Node::Join(process_subsplit(parent_subsplit),
                      process_subsplit(parent_subsplit.RotateSubsplit()));
  }

  void CheckSBNMapsAvailable() {
    if (!indexer_.size() || !index_to_child_.size() ||
        parent_to_range_.size() || !rootsplits_.size()) {
      std::cout << "Please call ProcessLoadedTrees to prepare your SBN maps.\n";
      abort();
    }
  }

  void SampleTrees(size_t count) {
    CheckSBNMapsAvailable();
    auto leaf_count = rootsplits_[0].size();
    auto edge_count = 2 * static_cast<int>(leaf_count) - 1;
    tree_collection_.trees_.clear();
    for (size_t i = 0; i < count; i++) {
      std::vector<double> branch_lengths(static_cast<size_t>(edge_count));
      tree_collection_.trees_.emplace_back(
          Tree(SampleTopology(), std::move(branch_lengths)));
    }
  }

  // This function returns a vector of vectors. Each vector corresponds to a
  // single virtual rooting of the tree, and gives the indices of sbn_probs_
  // corresponding to PCSSs that are present in the given topology.
  std::pair<SizeVector, SizeVectorVector> IndexerRepresentationOfTopology(
      Node::NodePtr topology) {
    auto tag_to_leafset = TagLeafSetMapOf(topology);
    auto leaf_count = topology->LeafCount();
    SizeVectorVector pcss_result(topology->Index());
    // First the hard part: the pcss_result.
    topology->PCSSPreOrder([
          &indexer_ = this->indexer_, &tag_to_leafset, &leaf_count, &pcss_result
    ](const Node *sister_node, bool sister_direction, const Node *focal_node,
      bool focal_direction,  //
      const Node *child0_node,
      bool child0_direction,  //
      const Node *child1_node, bool child1_direction,
      const Node *virtual_root_clade) {
      // Start by making the bitset representation of this PCSS.
      Bitset bitset(3 * leaf_count, false);
      bitset.CopyFrom(tag_to_leafset.at(sister_node->Tag()), 0,
                      sister_direction);
      bitset.CopyFrom(tag_to_leafset.at(focal_node->Tag()), leaf_count,
                      focal_direction);
      auto child0_bitset = tag_to_leafset.at(child0_node->Tag());
      if (child0_direction) child0_bitset.flip();
      auto child1_bitset = tag_to_leafset.at(child1_node->Tag());
      if (child1_direction) child1_bitset.flip();
      bitset.CopyFrom(std::min(child0_bitset, child1_bitset), 2 * leaf_count,
                      false);
      auto indexer_position = indexer_.at(bitset);
      const auto &focal_index = focal_node->Index();
      if (sister_node == focal_node) {
        // We are in the bidirectional edge situation.
        assert(focal_index < pcss_result.size());
        // Rooting at the present edge will indeed lead to the given PCSS.
        pcss_result[focal_index].push_back(indexer_position);
      } else {
        // The only time the virtual root clade should be nullptr should be when
        // sister_node == focal_node, but we check anyhow.
        assert(virtual_root_clade != nullptr);
        // Virtual-rooting on every edge in the sister will also lead to this
        // PCSS, because then the root will be "above" the PCSS.
        virtual_root_clade->ConditionalPreOrder(
            [&pcss_result, &indexer_position, &sister_node,
             &focal_node](const Node *node) {
              if (node == sister_node || node == focal_node) {
                return false;
              }
              pcss_result[node->Index()].push_back(indexer_position);
              return true;
            });
      }
    });
    // Second, the rootsplits.
    SizeVector rootsplit_result(topology->Index() + 1);
    topology->PreOrder([
      &rootsplit_result, &tag_to_leafset, &indexer = this->indexer_
    ](const Node *node) {
      Bitset rootsplit = tag_to_leafset.at(node->Tag());
      rootsplit.Minorize();
      rootsplit_result[node->Index()] = indexer.at(rootsplit);
    });
  return std::pair<SizeVector, SizeVectorVector>(rootsplit_result, pcss_result);
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
    StringUInt32PairMap parent_to_range_str;
    for (const auto &iter : parent_to_range_) {
      assert(parent_to_range_str.insert({iter.first.ToString(), iter.second})
                 .second);
    }
    assert(parent_to_range_str.insert({"rootsplit", {0, rootsplit_index_end_}})
               .second);
    return std::tie(indexer_str, parent_to_range_str);
  }

  // This function is really just for testing-- it recomputes from scratch.
  std::pair<StringUInt32Map, StringPCSSMap> SplitCounters() {
    auto counter = tree_collection_.TopologyCounter();
    return {StringUInt32MapOf(RootsplitCounterOf(counter)),
            StringPCSSMapOf(PCSSCounterOf(counter))};
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

  void CheckDataLoaded() {
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
  }

  void CheckBeagleDimensions() {
    CheckDataLoaded();
    if (beagle_instances_.size() == 0) {
      std::cerr << "Call MakeBeagleInstances to make some instances for "
                   "likelihood computation.\n";
      abort();
    }
    if (alignment_.SequenceCount() != beagle_leaf_count_ ||
        alignment_.Length() != beagle_site_count_) {
      std::cerr << "Alignment dimensions for current BEAGLE instances do not "
                   "match current alignment. Call MakeBeagleInstances again.\n";
      abort();
    }
  }

  void MakeBeagleInstances(int instance_count) {
    // Start by clearing out any existing instances.
    FinalizeBeagleInstances();
    CheckDataLoaded();
    beagle_leaf_count_ = alignment_.SequenceCount();
    beagle_site_count_ = alignment_.Length();
    SitePattern site_pattern(alignment_, tree_collection_.TagTaxonMap());
    for (auto i = 0; i < instance_count; i++) {
      auto beagle_instance = beagle::CreateInstance(site_pattern);
      beagle::SetJCModel(beagle_instance);
      beagle_instances_.push_back(beagle_instance);
      beagle::PrepareBeagleInstance(beagle_instance, tree_collection_,
                                    site_pattern);
    }
  }

  std::vector<double> LogLikelihoods() {
    CheckBeagleDimensions();
    return beagle::LogLikelihoods(beagle_instances_, tree_collection_);
  }

  std::vector<std::pair<double, std::vector<double>>> BranchGradients() {
    return beagle::BranchGradients(beagle_instances_, tree_collection_);
  }
};

// Here we initialize our static random number generator.
std::random_device SBNInstance::random_device_;
std::mt19937 SBNInstance::random_generator_(SBNInstance::random_device_());

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("libsbn") {
  SBNInstance inst("charlie");
  inst.ReadNewickFile("data/hello.nwk");
  inst.ReadFastaFile("data/hello.fasta");
  inst.MakeBeagleInstances(2);
  for (auto ll : inst.LogLikelihoods()) {
    CHECK_LT(abs(ll - -84.852358), 0.000001);
  }
  // Reading one file after another checks that we've cleared out state.
  inst.ReadNewickFile("data/five_taxon.nwk");
  inst.ProcessLoadedTrees();
  auto indexer_test_topology =
      // (2,(1,3)5,(0,4)6)7
      Node::OfParentIndexVector({6, 5, 7, 5, 6, 7, 7});
  auto representation =
      inst.IndexerRepresentationOfTopology(indexer_test_topology);
  std::cout << representation;
  //  SizeVectorVector correct_representation({{54, 31, 15},
  //                                           {21, 13, 22},
  //                                           {74, 42, 61},
  //                                           {21, 12, 63},
  //                                           {27, 31, 57},
  //                                           {21, 75, 76},
  //                                           {70, 73, 31}});
  //  CHECK_EQ(representation, correct_representation);

  inst.ReadNexusFile("data/DS1.subsampled_10.t");
  inst.ReadFastaFile("data/DS1.fasta");
  inst.MakeBeagleInstances(2);
  auto likelihoods = inst.LogLikelihoods();
  std::vector<double> pybeagle_likelihoods(
      {-14582.995273982739, -6911.294207416366, -6916.880235529542,
       -6904.016888831189, -6915.055570693576, -6915.50496696512,
       -6910.958836661867, -6909.02639968063, -6912.967861935749,
       -6910.7871105783515});
  for (size_t i = 0; i < likelihoods.size(); i++) {
    CHECK_LT(abs(likelihoods[i] - pybeagle_likelihoods[i]), 0.00011);
  }

  auto gradients = inst.BranchGradients();
  // Test the log likelihoods.
  for (size_t i = 0; i < likelihoods.size(); i++) {
    CHECK_LT(abs(gradients[i].first - pybeagle_likelihoods[i]), 0.00011);
    std::cout << gradients[i].first << " " << pybeagle_likelihoods[i]
              << std::endl;
  }
  // Test the gradients for the last tree.
  auto last = gradients.back();
  std::sort(last.second.begin(), last.second.end());
  // Zeros are for the root and one of the descendants of the root.
  std::vector<double> physher_gradients = {
      -904.18956, -607.70500, -562.36274, -553.63315, -542.26058, -539.64210,
      -463.36511, -445.32555, -414.27197, -412.84218, -399.15359, -342.68038,
      -306.23644, -277.05392, -258.73681, -175.07391, -171.59627, -168.57646,
      -150.57623, -145.38176, -115.15798, -94.86412,  -83.02880,  -80.09165,
      -69.00574,  -51.93337,  0.00000,    0.00000,    16.17497,   20.47784,
      58.06984,   131.18998,  137.10799,  225.73617,  233.92172,  253.49785,
      255.52967,  259.90378,  394.00504,  394.96619,  396.98933,  429.83873,
      450.71566,  462.75827,  471.57364,  472.83161,  514.59289,  650.72575,
      888.87834,  913.96566,  927.14730,  959.10746,  2296.55028};
  for (size_t i = 0; i < last.second.size(); i++) {
    CHECK_LT(abs(last.second[i] - physher_gradients[i]), 0.0001);
  }
}
#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_LIBSBN_HPP_
