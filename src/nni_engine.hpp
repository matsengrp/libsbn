// Copyright 2019-2022 bito project contributors.
// bito is free software under the GPLv3; see LICENSE file for details.
//
// NNI Engine
// This engine is used to explore the topological space surrounding the
// subsplitDAG using the smallest topological change of the NNI (Nearest Neighbor
// Interchange). The engine functions by finding all the NNIs adjacent to the DAG,
// evaluating them by their likelihoods, then choosing to accept or reject each NNI
// based on a filtering criteria.  For each NNI added to the DAG, each new edge added
// also adds new adjacent NNIs to the DAG.  This process is repeated until no NNIs pass
// the filtering criteria.

#pragma once

#include "gp_engine.hpp"
#include "tp_engine.hpp"
#include "gp_dag.hpp"

#include "bitset.hpp"
#include "subsplit_dag.hpp"
#include "nni_operation.hpp"
#include "graft_dag.hpp"
#include "sugar.hpp"
#include "gp_operation.hpp"
#include "reindexer.hpp"

using NNIDoubleMap = std::map<NNIOperation, double>;

class NNIEngine {
 public:
  // Constructors:
  NNIEngine(GPDAG &dag, std::optional<GPEngine *> gp_engine = std::nullopt,
            std::optional<TPEngine *> tp_engine = std::nullopt);

  // ** Access

  // Get Reference of DAG.
  const GPDAG &GetGPDAG() const { return dag_; };
  // Get Reference of GraftDAG.
  GraftDAG &GetGraftDAG() const { return *graft_dag_.get(); };
  // Get Reference of GPEngine.
  const GPEngine &GetGPEngine() const { return *gp_engine_; }
  bool &IsUsingGPEngine() { return gp_engine_ != nullptr; }
  // Get Reference of TPEngine.
  const TPEngine &GetTPEngine() const { return *tp_engine_; }
  bool &IsUsingTPEngine() { return tp_engine_ != nullptr; }
  // Get Adjacent NNIs to DAG.
  const NNISet &GetAdjacentNNIs() const { return adjacent_nnis_; };
  // Get Number of Adjacent NNIs.
  size_t GetAdjacentNNICount() const { return adjacent_nnis_.size(); };
  // Get NNIs that have been accepted into the DAG on current pass.
  const NNISet &GetAcceptedNNIs() const { return accepted_nnis_; };
  // Get all past NNIs that have been accepted into the DAG during run.
  const NNISet &GetPastAcceptedNNIs() const { return accepted_past_nnis_; };
  // Get NNIs that have been rejected from the DAG on current pass.
  const NNISet &GetRejectedNNIs() const { return rejected_nnis_; };
  // Get all past NNIs that have been rejected from the DAG during run.
  const NNISet &GetPastRejectedNNIs() const { return rejected_past_nnis_; };
  // Get Number of Accepted NNIs during the current pass.
  size_t GetAcceptedNNICount() const { return accepted_nnis_.size(); };
  // Get Map of proposed NNIs with their score.
  const NNIDoubleMap &GetScoredNNIs() const { return scored_nnis_; };
  // Get number of runs of NNI engine.
  size_t GetSweepCount() const { return sweep_count_; };

  // Set GP Engine.
  void SetGPEngine(GPEngine *gp_engine) { gp_engine_ = gp_engine; }
  // Set TP Engine.
  void SetTPEngine(TPEngine *tp_engine) { tp_engine_ = tp_engine; }

  // ** Runners
  // These start the engine, which procedurally ranks and adds (and maybe removes) NNIs
  // to the DAG, until some termination criteria has been satisfied.

  // Primary Runner for NNI Engine.
  void Run();
  // Initialization step run before loop.
  void RunInit();
  // Step that finds adjacent NNIs, evaluates, tehn accepts or rejects them.
  void RunMainLoop();
  // Step that runs at the end of each loop, preps for next loop.
  void RunPostLoop();

  // ** Filter Functions

  // Initialization step for filter before beginning
  using StaticFilterInitFunction =
      std::function<void(NNIEngine &, GPEngine &, TPEngine &, GraftDAG &)>;
  // Update step for any filter computation work to be done at start of each sweep.
  using StaticFilterUpdateFunction =
      std::function<void(NNIEngine &, GPEngine &, TPEngine &, GraftDAG &)>;
  // Function template for computational evaluation to be performed on an adjacent NNI.
  using StaticFilterEvaluateFunction = std::function<double(
      NNIEngine &, GPEngine &, TPEngine &, GraftDAG &, const NNIOperation &)>;
  // Function template for processing an adjacent NNI to be accepted or rejected.
  using StaticFilterProcessFunction = std::function<bool(
      NNIEngine &, GPEngine &, TPEngine &, GraftDAG &, const NNIOperation &)>;

  // Set to evaluate all NNIs to 0.
  void SetNoEvaluate();
  // Set filter to accept/deny all adjacent NNIs.
  void SetNoFilter(const bool set_nni_to_pass = true);
  // Set cutoff filter to constant cutoff.
  void SetScoreCutoff(const double score_cutoff);

  // ** Filtering

  // Initialize filter before first NNI sweep.
  void FilterInit();
  // Update filter parameters for each NNI sweep (before evaluation).
  void FilterPreUpdate();
  // Perform computation on each Adjacent NNI.
  void FilterEvaluateAdjacentNNIs();
  // Update filter parameters for each NNI sweep (after evaluation).
  void FilterPostUpdate();
  // Apply the filtering method to determine whether each Adjacent NNI will be added to
  // Accepted NNI or Rejected NNI.
  void FilterProcessAdjacentNNIs();

  // ** Filtering Schemes

  // Set custom filtering scheme functions.
  void SetFilteringScheme(
      std::optional<StaticFilterInitFunction> filter_init_fn = std::nullopt,
      std::optional<StaticFilterUpdateFunction> filter_pre_update_fn = std::nullopt,
      std::optional<StaticFilterEvaluateFunction> filter_eval_fn = std::nullopt,
      std::optional<StaticFilterUpdateFunction> filter_post_update_fn = std::nullopt,
      std::optional<StaticFilterProcessFunction> filter_process_fn = std::nullopt);

  // Set filtering scheme to use GP likelihoods.
  void SetGPLikelihoodFilteringScheme(const double score_cutoff);
  // Set filtering scheme to use TP likelihoods.
  void SetTPLikelihoodFilteringScheme(const double score_cutoff);
  // Set filtering scheme to use TP parsimony.
  void SetTPParsimonyFilteringScheme(const double score_cutoff);

  // ** Key Indexing
  // These function finds mapping from current NNIs contained in DAG to proposed NNI.

  // Map for storing important key indices for NNI Likelihood computation.
  enum class KeyIndex : size_t {
    Parent_Id,
    Child_Id,
    Edge,
    Parent_RHat,
    Parent_RFocal,
    Parent_PHatSister,
    Child_P,
    Child_PHatLeft,
    Child_PHatRight,
  };
  static const size_t KeyIndexCount = 9;
  using KeyIndexIterator =
      EnumIterator<KeyIndex, KeyIndex::Parent_Id, KeyIndex::Child_PHatRight>;
  using KeyIndexMap = EnumArray<KeyIndex, KeyIndexCount, size_t>;
  using KeyIndexMapPair = std::pair<KeyIndexMap, KeyIndexMap>;
  using KeyIndexPairArray = std::array<std::pair<KeyIndex, KeyIndex>, 3>;

  // Translate NNIClade type to corresponding PHat PLV from KeyIndex type.
  using NNIClade = NNIOperation::NNIClade;
  static KeyIndex NNICladeToPHatPLV(NNIClade clade_type);

  // Builds an array containing a mapping of plvs from Pre-NNI to Post-NNI, according to
  // remapping of the NNI clades (sister, right_child, left_child).
  static KeyIndexPairArray BuildKeyIndexTypePairsFromPreNNIToPostNNI(
      const NNIOperation &pre_nni, const NNIOperation &post_nni);

  // Create map of key indices for given NNIOperation needed for computing NNI
  // Likelihoods.
  template <typename DAGType>
  static KeyIndexMap BuildKeyIndexMapForNNI(const NNIOperation &nni, const DAGType &dag,
                                            const size_t node_count);
  KeyIndexMap BuildKeyIndexMapForNNI(const NNIOperation &nni,
                                     const size_t node_count) const;
  // Create map of key indices for Post-NNI, using Pre-NNI map as a reference.
  template <typename DAGType>
  static KeyIndexMap BuildKeyIndexMapForPostNNIViaReferencePreNNI(
      const NNIOperation &pre_nni, const NNIOperation &post_nni,
      const KeyIndexMap &pre_key_idx, const DAGType &dag);
  KeyIndexMap BuildKeyIndexMapForPostNNIViaReferencePreNNI(
      const NNIOperation &pre_nni, const NNIOperation &post_nni,
      const KeyIndexMap &pre_key_idx);

  // ** Scoring via GP Likelihood
  // Functions for computation or estimating likelihood of NNIs.

  // Fetches Pre-NNI data to prep Post-NNI for likelihood computation. Method stores
  // intermediate values in the GPEngine temp space (expects GPEngine has already been
  // resized).
  KeyIndexMapPair PassGPEngineDataFromPreNNIToPostNNIViaReference(
      const NNIOperation &pre_nni, const NNIOperation &post_nni,
      const size_t nni_count = 0, const bool use_unique_temps = false);
  // Fetches Data from Pre-NNI for Post-NNI. Can be used for initial values when moving
  // accepted NNIs from Graft to Host DAG.
  KeyIndexMapPair PassGPEngineDataFromPreNNIToPostNNIViaCopy(
      const NNIOperation &pre_nni, const NNIOperation &post_nni);

  // Build GPOperation vector for computing NNI Likelihood. For NNIs that are not in
  // the DAG, must create KeyIndexMap through
  // PassGPEngineDataFromPreNNIToPostNNIViaReference.
  GPOperationVector BuildGPOperationsForNNILikelihood(
      const NNIOperation &nni, const KeyIndexMap &nni_key_idx) const;
  // Build GPOperation vector for computing NNI Likelihood for all adjacent NNIs.
  GPOperationVector BuildGPOperationsForAdjacentNNILikelihoods(
      const bool via_reference = true);
  // Grow engine to handle computing NNI Likelihoods for all adjacent NNIs.
  // Option to grow engine for computing likelihoods via reference or via copy. If
  // computing via reference, option whether to use unique temporaries (for testing and
  // computing in parallel).
  void GrowGPEngineForAdjacentNNILikelihoods(const bool via_reference = true,
                                             const bool use_unique_temps = false);
  // Performs entire GP likelihood computation for all Adjacent NNIs.
  // Allocates necessary extra space on GPEngine, generates and processes
  // GPOperationVector, then retrieves results and stores in Scored NNIs.
  void ScoreAdjacentNNIsByGPLikelihood();

  // ** Scoring via TP Likelihood

  // Performs entire TP likelihood computation for all Adjacent NNIs.
  // Allocates necessary extra space on TPEngine, computes and retrieves results and
  // stores in Scored NNIs.
  void ScoreAdjacentNNIsByTPLikelihood();

  // ** Scoring via TP Parsimony

  void ScoreAdjacentNNIsByTPParsimony();

  // ** DAG & GPEngine Maintenance

  // Initial GPEngine for use with GraftDAG.
  void InitGPEngine();
  // Populate PLVs so that
  void PrepGPEngineForLikelihoods();
  // Add all Accepted NNIs to Main DAG.
  void AddAcceptedNNIsToDAG();
  // Add all Adjacent NNIs to Graft DAG.
  void GraftAdjacentNNIsToDAG();
  // Remove all NNIs from Graft DAG.
  void RemoveAllGraftedNNIsFromDAG();

  // ** NNI Maintenance
  // These maintain NNIs to stay consistent with the state of associated GraftDAG.

  // Freshly synchonizes NNISet to match the current state of its DAG. Wipes old NNI
  // data and finds all all parent/child pairs adjacent to DAG by iterating over all
  // internal edges in the DAG. (For each internal edges, two NNIs are possible.)
  void SyncAdjacentNNIsWithDAG();
  // Updates NNI Set after given parent/child node pair have been added to the DAG.
  // Removes pair from NNI Set and adds adjacent pairs coming from newly created edges.
  void UpdateAdjacentNNIsAfterDAGAddNodePair(const NNIOperation &nni);
  void UpdateAdjacentNNIsAfterDAGAddNodePair(const Bitset &parent_bitset,
                                             const Bitset &child_bitset);
  // Adds all NNIs from all (node_id, other_id) pairs, where other_id's are elements of
  // the adjacent_node_ids vector. is_edge_leafward tells whether node_id is the child
  // or parent. is_edge_edgeclade determines which side of parent the child descends
  // from.
  void AddAllNNIsFromNodeVectorToAdjacentNNIs(const NodeId node_id,
                                              const SizeVector &adjacent_node_ids,
                                              const bool is_edge_on_left,
                                              const bool is_edge_leafward);
  // Based on given input NNIOperation, produces the two possible output NNIOperations
  // and adds those results to the NNI Set (if results are not a member of the DAG or
  // NNI Set).
  void SafeAddOutputNNIsToAdjacentNNIs(const Bitset &parent_bitset,
                                       const Bitset &child_bitset,
                                       const bool is_edge_on_left);

  // Add score to given NNI.
  void AddScoreForNNI(const NNIOperation &nni, const double score);
  // Get score for given NNI.
  double GetScoreForNNI(const NNIOperation &nni) const;

  // Update adjacent NNIs at end of current sweep. If not re-evaluating rejected NNIs,
  // then current adjacent NNIs are removed. Then NNIs that are adjacent to last sweep's
  // accepted NNIs are added.
  void UpdateAdjacentNNIs(const bool reevaluate_rejected_nnis = false);
  // Remove all accepted NNIs and optionally save to past NNIs.
  void UpdateAcceptedNNIs(const bool save_past_nnis = true);
  // Remove all rejected NNIs and optionally save to past NNIs.
  void UpdateRejectedNNIs(const bool save_past_nnis = true);
  // Remove all scored NNIs and optionally save to past NNIs.
  void UpdateScoredNNIs(const bool save_past_nnis = false);
  // Reset all NNIs, current and past.
  void ResetAllNNIs();

  // ** Miscellaneous

  NNIOperation FindNNINeighborInDAG(const NNIOperation &nni);

 private:
  // ** Access

  // Get Reference of GPEngine.
  GPEngine &GetGPEngine() { return *gp_engine_; }
  // Get Reference of GPEngine.
  TPEngine &GetTPEngine() { return *tp_engine_; }

 private:
  // Un-owned reference DAG.
  GPDAG &dag_;
  // For adding temporary NNIs to DAG.
  std::unique_ptr<GraftDAG> graft_dag_;
  // Un-owned reference GPEngine.
  GPEngine *gp_engine_ = nullptr;
  // Un-owned reference TPEngine.
  TPEngine *tp_engine_ = nullptr;
  // Tracks modifications to the DAG.
  Reindexer node_reindexer_;
  Reindexer edge_reindexer_;

  // Set of NNIs to be evaluated, which are a single NNI.
  NNISet adjacent_nnis_;
  // NNIs which have passed the filtering threshold, to be added to the DAG.
  NNISet accepted_nnis_;
  // NNIs which have been accepted in a previous sweep of the search.
  NNISet accepted_past_nnis_;
  // NNIs which have failed the filtering threshold, not to be added to the DAG.
  NNISet rejected_nnis_;
  // NNIs which have been rejected in a previous sweep of the search.
  NNISet rejected_past_nnis_;

  // Map of adjacent NNIs to their score.
  NNIDoubleMap scored_nnis_;
  // Map of previous rejected NNIs to their score.
  NNIDoubleMap scored_past_nnis_;

  // Steps of filtering scheme.
  StaticFilterInitFunction filter_init_fn_ = nullptr;
  StaticFilterUpdateFunction filter_pre_update_fn_ = nullptr;
  StaticFilterEvaluateFunction filter_eval_fn_ = nullptr;
  StaticFilterUpdateFunction filter_post_update_fn_ = nullptr;
  StaticFilterProcessFunction filter_process_fn_ = nullptr;

  // Count number of loops executed by engine.
  size_t sweep_count_ = 0;
};
