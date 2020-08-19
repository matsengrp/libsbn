// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_ROOTED_SBN_INSTANCE_HPP_
#define SRC_ROOTED_SBN_INSTANCE_HPP_

#include "generic_sbn_instance.hpp"
#include "rooted_sbn_support.hpp"

using PreRootedSBNInstance = GenericSBNInstance<RootedTreeCollection, RootedSBNSupport>;
template class GenericSBNInstance<RootedTreeCollection, RootedSBNSupport>;

class RootedSBNInstance : public PreRootedSBNInstance {
 public:
  using PreRootedSBNInstance::PreRootedSBNInstance;

  // ** SBN-related items

  // Turn an IndexerRepresentation into a string representation of the underying
  // bitsets. This is really just so that we can make a test of indexer
  // representations.
  StringSet StringIndexerRepresentationOf(
      const RootedIndexerRepresentation& indexer_representation) const;
  StringSet StringIndexerRepresentationOf(const Node::NodePtr& topology,
                                          size_t out_of_sample_index) const;

  // ** Phylogenetic likelihood

  std::vector<double> LogLikelihoods();
  // For each loaded tree, return the phylogenetic gradient.
  std::vector<RootedPhyloGradient> PhyloGradients();

  // ** I/O

  void ReadNewickFile(std::string fname);
  void ReadNexusFile(std::string fname);
};

#ifdef DOCTEST_LIBRARY_INCLUDED

#include "doctest_constants.hpp"

// Centered finite difference approximation of the derivative wrt rate.
std::vector<double> derivative_strict_clock(RootedSBNInstance& inst) {
  double eps = 0.00000001;
  std::vector<double> rates;
  std::vector<double> gradients;

  for (auto& tree : inst.tree_collection_.trees_) {
    rates.push_back(tree.rates_[0]);
    tree.rates_.assign(tree.rates_.size(), rates.back() - eps);
  }
  auto lm = inst.LogLikelihoods();

  int i = 0;
  for (auto& tree : inst.tree_collection_.trees_) {
    tree.rates_.assign(tree.rates_.size(), rates[i++] + eps);
  }
  auto lp = inst.LogLikelihoods();

  for (size_t index = 0; index < lm.size(); index++) {
    gradients.push_back((lp[index] - lm[index]) / (2. * eps));
  }
  return gradients;
}

// Centered finite difference approximation of the derivative wrt to each rate.
std::vector<std::vector<double>> derivative_relaxed_clock(RootedSBNInstance& inst) {
  double eps = 0.00000001;
  std::vector<std::vector<double>> gradients;
  std::vector<double> lp;
  std::vector<double> lm;
  size_t edge_count = inst.TaxonCount() * 2 - 2;

  for (size_t index = 0; index < edge_count; index++) {
    std::vector<double> gradient;
    std::vector<double> rates;
    for (int i = 0; i < inst.tree_collection_.TreeCount(); i++) {
      double value = inst.tree_collection_.trees_[i].rates_[index];
      rates.push_back(value);
      inst.tree_collection_.trees_[i].rates_[index] = rates.back() - eps;
    }
    lm = inst.LogLikelihoods();

    for (size_t i = 0; i < inst.tree_collection_.TreeCount(); i++) {
      inst.tree_collection_.trees_[i].rates_[index] = rates[i] + eps;
    }
    lp = inst.LogLikelihoods();

    for (size_t i = 0; i < inst.tree_collection_.TreeCount(); i++) {
      inst.tree_collection_.trees_[i].rates_[index] = rates[i];
      gradient.push_back((lp[i] - lm[i]) / (2. * eps));
    }

    gradients.push_back(gradient);
  }
  return gradients;
}

TEST_CASE("RootedSBNInstance: subsplit support and TrainSimpleAverage") {
  RootedSBNInstance inst("charlie");
  inst.ReadNewickFile("data/five_taxon_rooted.nwk");
  inst.ProcessLoadedTrees();
  auto pretty_indexer = inst.PrettyIndexer();
  StringSet pretty_indexer_set{pretty_indexer.begin(), pretty_indexer.end()};
  // The indexer_ is to index the sbn_parameters_. Note that neither of these
  // data structures attempt to catalog the complete collection of rootsplits or
  // PCSPs, but just those that are present in the the input trees.
  //
  // The indexer_ and sbn_parameters_ are laid out as follows (I'll just call it
  // the "index" in what follows). Say there are rootsplit_count rootsplits in
  // the support.
  // The first rootsplit_count entries of the index are assigned to the
  // rootsplits (again, those rootsplits that are present for some rooting of
  // the unrooted input trees). The rest of the entries of the index are laid out as
  // blocks of parameters for PCSPs that share the same parent. Take a look at the
  // description of PCSP bitsets (and the unit tests) in bitset.hpp to understand the
  // notation used here.
  //
  // In contrast to the unrooted case, we can write out the pretty indexer here and
  // verify it by hand. There is the block structure in which the two children of
  // 10000|01111 are grouped together.
  StringSet correct_pretty_indexer_set{
      "00111",              // ((x0,x1),(x2,(x3,x4)))
      "01111",              // (x0,(((x1,x3),x2),x4)) and ((x1,((x2,x4),x3)),x0)
      "00010",              // (x3,((x0,(x4,x1)),x2))
      "00100|01010|00010",  // ((x1,x3),x2)
      "00111|11000|01000",  // ((x0,x1),(x2,(x3,x4)))
      "00100|00011|00001",  // (x2,(x3,x4))
      "11000|00111|00011",  // ((x0,x1),(x2,(x3,x4)))
      "00100|11001|01001",  // ((x0,(x4,x1)),x2)
      "10000|01001|00001",  // (x0,(x4,x1))
      "01000|00111|00010",  // (x1,((x2,x4),x3))
      "10000|01111|00001",  // (x0,(((x1,x3),x2),x4))
      "10000|01111|00111",  // ((x1,((x2,x4),x3)),x0)
      "00010|00101|00001",  // ((x2,x4),x3)
      "00001|01110|00100",  // (((x1,x3),x2),x4)
      "00010|11101|00100"   // (x3,((x0,(x4,x1)),x2))
  };
  CHECK_EQ(pretty_indexer_set, correct_pretty_indexer_set);

  // Test of rooted IndexerRepresentationOf.
  // Topology is ((0,1),(2,(3,4)));, or with internal nodes ((0,1)5,(2,(3,4)6)7)8;
  auto indexer_test_rooted_topology = Node::OfParentIdVector({5, 5, 7, 6, 6, 8, 7, 8});
  auto correct_rooted_indexer_representation = StringSet(
      {"00111", "11000|00111|00011", "00100|00011|00001", "00111|11000|01000"});
  CHECK_EQ(inst.StringIndexerRepresentationOf(indexer_test_rooted_topology,
                                              out_of_sample_index),
           correct_rooted_indexer_representation);

  inst.TrainSimpleAverage();
  StringVector correct_taxon_names({"x0", "x1", "x2", "x3", "x4"});
  CHECK_EQ(inst.SBNSupport().TaxonNames(), correct_taxon_names);
  StringDoubleVector correct_parameters({{"00111", 0.25},
                                         {"01111", 0.5},
                                         {"00010", 0.25},
                                         {"00100|01010|00010", 1},
                                         {"00111|11000|01000", 1},
                                         {"00100|00011|00001", 1},
                                         {"11000|00111|00011", 1},
                                         {"00100|11001|01001", 1},
                                         {"10000|01001|00001", 1},
                                         {"01000|00111|00010", 1},
                                         {"10000|01111|00001", 0.5},
                                         {"10000|01111|00111", 0.5},
                                         {"00010|00101|00001", 1},
                                         {"00001|01110|00100", 1},
                                         {"00010|11101|00100", 1}});
  CHECK_EQ(correct_parameters, inst.PrettyIndexedSBNParameters());
}

TEST_CASE("RootedSBNInstance: gradients") {
  RootedSBNInstance inst("charlie");
  inst.ReadNewickFile("data/fluA.tree");
  inst.ReadFastaFile("data/fluA.fa");
  PhyloModelSpecification simple_specification{"JC69", "constant", "strict"};
  inst.PrepareForPhyloLikelihood(simple_specification, 1);
  for (auto& tree : inst.tree_collection_.trees_) {
    tree.rates_.assign(tree.rates_.size(), 0.001);
  }

  auto likelihood = inst.LogLikelihoods();
  double physher_ll = -4777.616349;
  CHECK_LT(fabs(likelihood[0] - physher_ll), 0.0001);

  auto gradients = inst.PhyloGradients();
  std::vector<double> physher_gradients = {
      -0.593654, 6.441290,   11.202945, 5.173924,  -0.904631, 2.731402,   3.157131,
      7.082914,  10.305417,  13.988206, 20.709336, 48.897993, 99.164949,  130.205747,
      17.314019, 21.033290,  -1.336335, 12.259822, 22.887291, 27.176564,  47.487426,
      3.637276,  12.955169,  15.315953, 83.254605, -3.806996, 105.385095, 4.874023,
      22.754466, 6.036534,   25.651478, 29.535185, 29.598789, 1.817247,   10.598685,
      76.259248, 56.481423,  10.679778, 6.587179,  3.330556,  -4.622247,  33.417304,
      63.415767, 188.809515, 23.540875, 17.421076, 1.222568,  22.372012,  34.239511,
      3.486115,  4.098873,   13.200954, 19.726890, 96.808738, 4.240029,   7.414585,
      48.871694, 3.488516,   82.969065, 9.009334,  8.032474,  3.981016,   6.543650,
      53.702423, 37.835952,  2.840831,  7.517186,  19.936861};
  for (size_t i = 0; i < physher_gradients.size(); i++) {
    CHECK_LT(fabs(gradients[0].ratios_root_height_[i] - physher_gradients[i]), 0.0001);
  }
  CHECK_LT(fabs(gradients[0].log_likelihood_ - physher_ll), 0.0001);
}

TEST_CASE("RootedSBNInstance: clock gradients") {
  RootedSBNInstance inst("charlie");
  inst.ReadNewickFile("data/fluA.tree");
  inst.ReadFastaFile("data/fluA.fa");
  PhyloModelSpecification simple_specification{"JC69", "constant", "strict"};
  inst.PrepareForPhyloLikelihood(simple_specification, 1);

  for (auto& tree : inst.tree_collection_.trees_) {
    tree.rates_.assign(tree.rates_.size(), 0.001);
  }

  auto likelihood = inst.LogLikelihoods();
  double physher_ll = -4777.616349;
  CHECK_LT(fabs(likelihood[0] - physher_ll), 0.0001);

  // Gradient with a strict clock.
  auto gradients_strict = inst.PhyloGradients();
  std::vector<double> gradients_strict_approx = derivative_strict_clock(inst);
  CHECK_LT(fabs(gradients_strict[0].clock_model_[0] - gradients_strict_approx[0]),
           0.001);
  CHECK_LT(fabs(gradients_strict[0].log_likelihood_ - physher_ll), 0.001);

  // Gradient with a "relaxed" clock.
  auto& tree = inst.tree_collection_.trees_[0];
  // Make a clock with some rate variation.
  for (size_t i = 0; i < tree.rates_.size(); i++) {
    tree.rates_[i] *= i % 3 + 1.0;
  }
  tree.rate_count_ = tree.rates_.size();

  auto gradients_relaxed = inst.PhyloGradients();
  auto gradients_relaxed_approx = derivative_relaxed_clock(inst);

  for (size_t j = 0; j < gradients_relaxed_approx.size(); j++) {
    CHECK_LT(
        fabs(gradients_relaxed[0].clock_model_[j] - gradients_relaxed_approx[j][0]),
        0.001);
  }
}

TEST_CASE("RootedSBNInstance: Weibull gradients") {
  RootedSBNInstance inst("charlie");
  inst.ReadNewickFile("data/fluA.tree");
  inst.ReadFastaFile("data/fluA.fa");
  PhyloModelSpecification simple_specification{"JC69", "weibull+4", "strict"};
  inst.PrepareForPhyloLikelihood(simple_specification, 1);

  for (auto& tree : inst.tree_collection_.trees_) {
    tree.rates_.assign(tree.rates_.size(), 0.001);
  }
  auto param_block_map = inst.GetPhyloModelParamBlockMap();
  param_block_map.at(WeibullSiteModel::shape_key_).setConstant(0.1);

  auto likelihood = inst.LogLikelihoods();
  double physher_ll = -4618.2062529058;
  CHECK_LT(fabs(likelihood[0] - physher_ll), 0.0001);

  // Gradient wrt Weibull site model.
  auto gradients = inst.PhyloGradients();
  double physher_gradient = -5.231329;
  CHECK_LT(fabs(gradients[0].site_model_[0] - physher_gradient), 0.001);
  CHECK_LT(fabs(gradients[0].log_likelihood_ - physher_ll), 0.001);
}

TEST_CASE("RootedSBNInstance: parsing dates") {
  RootedSBNInstance inst("charlie");
  inst.ReadNexusFile("data/test_beast_tree_parsing.nexus");
  std::vector<double> dates;
  for (auto [tag, date] : inst.tree_collection_.tag_date_map_) {
    dates.push_back(date);
  }
  std::sort(dates.begin(), dates.end());
  CHECK_EQ(dates[0], 0);
  CHECK_EQ(dates.back(), 80.0);
}

#endif  // DOCTEST_LIBRARY_INCLUDED
#endif  // SRC_ROOTED_SBN_INSTANCE_HPP_
