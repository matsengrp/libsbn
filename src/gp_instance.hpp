// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_GP_INSTANCE_HPP_
#define SRC_GP_INSTANCE_HPP_

#include "gp_dag.hpp"
#include "gp_engine.hpp"
#include "rooted_tree_collection.hpp"
#include "site_pattern.hpp"

class GPInstance {
 public:
  explicit GPInstance(const std::string &mmap_file_path)
      : mmap_file_path_(mmap_file_path) {
    if (mmap_file_path.empty()) {
      Failwith("GPInstance needs a legal path as a constructor argument.");
    }
  };
  void PrintStatus();

  void ReadFastaFile(const std::string &fname);
  void ReadNewickFile(const std::string &fname);
  void ReadNexusFile(const std::string &fname);

  void MakeEngine(double rescaling_threshold = GPEngine::default_rescaling_threshold_);
  GPEngine *GetEngine() const;
  bool HasEngine() const;
  void PrintDAG();
  void PrintGPCSPIndexer();
  void ProcessOperations(const GPOperationVector &operations);
  void HotStartBranchLengths();
  void EstimateSBNParameters();
  void EstimateBranchLengths(double tol, size_t max_iter);
  void ResetMarginalLikelihoodAndPopulatePLVs();
  void ComputeLikelihoods();
  RootedTreeCollection GenerateCompleteRootedTreeCollection();

  // #273: A lot of code duplication here with things in SBNInstance.
  StringVector PrettyIndexer() const;
  StringDoubleVector PrettyIndexedSBNParameters();
  void SBNParametersToCSV(const std::string &file_path);

 private:
  std::string mmap_file_path_;
  Alignment alignment_;
  std::unique_ptr<GPEngine> engine_;
  RootedTreeCollection tree_collection_;
  GPDAG dag_;
  static constexpr size_t plv_count_per_node_ = 6;

  // A vector that contains all of the SBN-related probabilities.
  EigenVectorXd sbn_parameters_;

  void ClearTreeCollectionAssociatedState();
  void CheckSequencesAndTreesLoaded() const;
  void ProcessLoadedTrees();

  void InitializeGPEngine();

  size_t GetGPCSPIndexForLeafNode(const Bitset &parent_subsplit, const Node *leaf_node);
};

#endif  // SRC_GP_INSTANCE_HPP_
