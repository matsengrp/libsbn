// Copyright 2019-2020 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "gp_engine.hpp"

#include "optimization.hpp"

GPEngine::GPEngine(SitePattern site_pattern, size_t plv_count, size_t gpcsp_count,
                   const std::string& mmap_file_path, double rescaling_threshold)
    : site_pattern_(std::move(site_pattern)),
      plv_count_(plv_count),
      rescaling_threshold_(rescaling_threshold),
      log_rescaling_threshold_(log(rescaling_threshold)),
      mmapped_master_plv_(mmap_file_path, plv_count_ * site_pattern_.PatternCount()) {
  Assert(plv_count_ > 0, "Zero PLV count in constructor of GPEngine.");
  plvs_ = mmapped_master_plv_.Subdivide(plv_count_);
  Assert(plvs_.size() == plv_count_,
         "Didn't get the right number of PLVs out of Subdivide.");
  Assert(plvs_.back().rows() == MmappedNucleotidePLV::base_count_ &&
             plvs_.back().cols() == site_pattern_.PatternCount(),
         "Didn't get the right shape of PLVs out of Subdivide.");
  plv_rescaling_counts_ = EigenVectorXi::Constant(plv_count_, 0);
  branch_lengths_.resize(gpcsp_count);
  branch_lengths_.setOnes();
  log_likelihoods_.resize(gpcsp_count);
  q_.resize(gpcsp_count);

  auto weights = site_pattern_.GetWeights();
  site_pattern_weights_ =
      Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(weights.data(), weights.size());

  InitializePLVsWithSitePatterns();
}

void GPEngine::operator()(const GPOperations::Zero& op) {
  plvs_.at(op.dest_).setZero();
  plv_rescaling_counts_(op.dest_) = 0;
}

void GPEngine::operator()(const GPOperations::SetToStationaryDistribution& op) {
  auto& plv = plvs_.at(op.dest_);
  for (size_t row_idx = 0; row_idx < plv.rows(); ++row_idx) {
    plv.row(row_idx).array() = stationary_distribution_(row_idx);
  }
  plv_rescaling_counts_(op.dest_) = 0;
}

void GPEngine::operator()(const GPOperations::IncrementWithWeightedEvolvedPLV& op) {
  SetTransitionMatrixToHaveBranchLength(branch_lengths_(op.gpcsp_));
  // We assume that we've done a PrepForMarginalization operation, and thus the
  // rescaling count for op.dest_ is the minimum of the rescaling counts among the
  // op.src_s. Thus this should be zero or negative:
  const int rescaling_difference =
      plv_rescaling_counts_(op.dest_) - plv_rescaling_counts_(op.src_);
  Assert(rescaling_difference <= 0,
         "Too much dest_ rescaling in IncrementWithWeightedEvolvedPLV");
  const double rescaling_factor =
      rescaling_difference == 0
          ? 1.
          : pow(rescaling_threshold_, static_cast<double>(rescaling_difference));
  // We are going to have evidence of reduced-precision arithmetic here because we are
  // adding together things of radically different rescaling amounts. This appears
  // unavoidable without special-purpose truncation code, which doesn't seem worthwhile.
  plvs_.at(op.dest_) +=
      rescaling_factor * q_(op.gpcsp_) * transition_matrix_ * plvs_.at(op.src_);
}

void GPEngine::operator()(const GPOperations::IncrementMarginalLikelihood& op) {
  Assert(plv_rescaling_counts_(op.stationary_) == 0,
         "Surprise! Rescaled stationary distribution in IncrementMarginalLikelihood");
  per_pattern_log_likelihoods_ =
      (plvs_.at(op.stationary_).transpose() * plvs_.at(op.p_))
          .diagonal()
          .array()
          .log() +
      LogRescalingFor(op.p_);

  log_likelihoods_[op.rootsplit_] =
      log(q_(op.rootsplit_)) + per_pattern_log_likelihoods_.dot(site_pattern_weights_);

  log_marginal_likelihood_ =
      NumericalUtils::LogAdd(log_marginal_likelihood_, log_likelihoods_[op.rootsplit_]);
}

void GPEngine::operator()(const GPOperations::Multiply& op) {
  plvs_.at(op.dest_).array() = plvs_.at(op.src1_).array() * plvs_.at(op.src2_).array();
  plv_rescaling_counts_(op.dest_) =
      plv_rescaling_counts_(op.src1_) + plv_rescaling_counts_(op.src2_);
  RescalePLVIfNeeded(op.dest_);
}

void GPEngine::operator()(const GPOperations::Likelihood& op) {
  SetTransitionMatrixToHaveBranchLength(branch_lengths_(op.dest_));
  PreparePerPatternLogLikelihoods(op.parent_, op.child_);
  log_likelihoods_[op.dest_] =
      log(q_[op.dest_]) + per_pattern_log_likelihoods_.dot(site_pattern_weights_);
}

void GPEngine::operator()(const GPOperations::OptimizeBranchLength& op) {
  BrentOptimization(op);
}

void GPEngine::operator()(const GPOperations::UpdateSBNProbabilities& op) {
  const size_t range_length = op.stop_ - op.start_;
  if (range_length == 1) {
    q_(op.start_) = 1.;
    return;
  }
  // else
  Eigen::VectorBlock<EigenVectorXd> segment =
      log_likelihoods_.segment(op.start_, range_length);
  const double log_norm = NumericalUtils::LogSum(segment);
  segment.array() -= log_norm;
  q_.segment(op.start_, range_length) = segment.array().exp();
}

void GPEngine::operator()(const GPOperations::PrepForMarginalization& op) {
  const size_t src_count = op.src_vector_.size();
  Assert(src_count > 0, "Empty src_vector in PrepForMarginalization");
  SizeVector src_rescaling_counts(src_count);
  for (size_t idx = 0; idx < src_count; ++idx) {
    src_rescaling_counts[idx] = plv_rescaling_counts_(op.src_vector_[idx]);
  }
  const auto min_rescaling_count =
      *std::min_element(src_rescaling_counts.begin(), src_rescaling_counts.end());
  plvs_.at(op.dest_).setZero();
  plv_rescaling_counts_(op.dest_) = min_rescaling_count;
}

void GPEngine::ProcessOperations(GPOperationVector operations) {
  for (const auto& operation : operations) {
    std::visit(*this, operation);
  }
}

void GPEngine::SetTransitionMatrixToHaveBranchLength(double branch_length) {
  diagonal_matrix_.diagonal() = (branch_length * eigenvalues_).array().exp();
  transition_matrix_ = eigenmatrix_ * diagonal_matrix_ * inverse_eigenmatrix_;
}

void GPEngine::SetTransitionAndDerivativeMatricesToHaveBranchLength(
    double branch_length) {
  diagonal_vector_ = (branch_length * eigenvalues_).array().exp();
  diagonal_matrix_.diagonal() = diagonal_vector_;
  transition_matrix_ = eigenmatrix_ * diagonal_matrix_ * inverse_eigenmatrix_;
  diagonal_matrix_.diagonal() = eigenvalues_.array() * diagonal_vector_.array();
  derivative_matrix_ = eigenmatrix_ * diagonal_matrix_ * inverse_eigenmatrix_;
}

void GPEngine::SetTransitionMatrixToHaveBranchLengthAndTranspose(double branch_length) {
  diagonal_matrix_.diagonal() = (branch_length * eigenvalues_).array().exp();
  transition_matrix_ =
      inverse_eigenmatrix_.transpose() * diagonal_matrix_ * eigenmatrix_.transpose();
}

void GPEngine::PrintPLV(size_t plv_idx) {
  for (const auto& row : plvs_[plv_idx].rowwise()) {
    std::cout << row << std::endl;
  }
  std::cout << std::endl;
}

DoublePair GPEngine::LogLikelihoodAndDerivative(
    const GPOperations::OptimizeBranchLength& op) {
  SetTransitionAndDerivativeMatricesToHaveBranchLength(branch_lengths_(op.gpcsp_));
  PreparePerPatternLogLikelihoods(op.rootward_, op.leafward_);
  // The prior is expressed using the current value of q_.
  // The phylogenetic component of the likelihood is weighted with the number of times
  // we see the site patterns.
  const double log_likelihood =
      log(q_(op.gpcsp_)) + per_pattern_log_likelihoods_.dot(site_pattern_weights_);

  // The per-site likelihood derivative is calculated in the same way as the per-site
  // likelihood, but using the derivative matrix instead of the transition matrix.
  // We first prepare two useful vectors _without_ likelihood rescaling, because the
  // rescalings cancel out in the ratio below.
  PrepareUnrescaledPerPatternLikelihoodDerivatives(op.rootward_, op.leafward_);
  PrepareUnrescaledPerPatternLikelihoods(op.rootward_, op.leafward_);
  // If l_i is the per-site likelihood, the derivative of log(l_i) is the derivative
  // of l_i divided by l_i.
  per_pattern_likelihood_derivative_ratios_ =
      per_pattern_likelihood_derivatives_.array() / per_pattern_likelihoods_.array();
  const double log_likelihood_derivative =
      per_pattern_likelihood_derivative_ratios_.dot(site_pattern_weights_);
  return {log_likelihood, log_likelihood_derivative};
}

void GPEngine::InitializePLVsWithSitePatterns() {
  for (auto& plv : plvs_) {
    plv.setZero();
  }
  size_t taxon_idx = 0;
  for (const auto& pattern : site_pattern_.GetPatterns()) {
    size_t site_idx = 0;
    for (const int symbol : pattern) {
      Assert(symbol >= 0, "Negative symbol!");
      if (symbol == MmappedNucleotidePLV::base_count_) {  // Gap character.
        plvs_.at(taxon_idx).col(site_idx).setConstant(1.);
      } else if (symbol < MmappedNucleotidePLV::base_count_) {
        plvs_.at(taxon_idx)(symbol, site_idx) = 1.;
      }
      site_idx++;
    }
    taxon_idx++;
  }
}

void GPEngine::RescalePLV(size_t plv_idx, int rescaling_count) {
  Assert(rescaling_count >= 0, "Negative rescaling count in RescalePLV.");
  if (rescaling_count == 0) {
    return;
  }
  // else
  plvs_.at(plv_idx) /= pow(rescaling_threshold_, static_cast<double>(rescaling_count));
  plv_rescaling_counts_(plv_idx) += rescaling_count;
}

void GPEngine::RescalePLVIfNeeded(size_t plv_idx) {
  double min_entry = plvs_.at(plv_idx).minCoeff();
  Assert(min_entry >= 0., "PLV with negative entry passed to RescalePLVIfNeeded");
  if (min_entry == 0) {
    return;
  }
  // else
  int rescaling_count = 0;
  while (min_entry < rescaling_threshold_) {
    min_entry /= rescaling_threshold_;
    rescaling_count++;
  }
  RescalePLV(plv_idx, rescaling_count);
}

double GPEngine::LogRescalingFor(size_t plv_idx) {
  return static_cast<double>(plv_rescaling_counts_(plv_idx)) * log_rescaling_threshold_;
}

void GPEngine::BrentOptimization(const GPOperations::OptimizeBranchLength& op) {
  auto negative_log_likelihood = [this, &op](double branch_length) {
    SetTransitionMatrixToHaveBranchLength(branch_length);
    PreparePerPatternLogLikelihoods(op.rootward_, op.leafward_);
    return -(log(q_[op.gpcsp_]) +
             per_pattern_log_likelihoods_.dot(site_pattern_weights_));
  };
  double current_branch_length = branch_lengths_(op.gpcsp_);
  double current_value = negative_log_likelihood(current_branch_length);
  const auto [branch_length, neg_log_likelihood] = Optimization::BrentMinimize(
      negative_log_likelihood, min_branch_length_, max_branch_length_,
      significant_digits_for_optimization_, max_iter_for_optimization_);

  // Numerical optimization sometimes yields new nllk > current nllk.
  // In this case, we reset the branch length to the previous value.
  if (neg_log_likelihood > current_value) {
    branch_lengths_(op.gpcsp_) = current_branch_length;
  } else {
    branch_lengths_(op.gpcsp_) = branch_length;
  }
}

void GPEngine::GradientAscentOptimization(
    const GPOperations::OptimizeBranchLength& op) {
  auto log_likelihood_and_derivative = [this, &op](double branch_length) {
    branch_lengths_(op.gpcsp_) = branch_length;
    return this->LogLikelihoodAndDerivative(op);
  };
  const auto [branch_length, log_likelihood] = Optimization::GradientAscent(
      log_likelihood_and_derivative, branch_lengths_(op.gpcsp_),
      relative_tolerance_for_optimization_, step_size_for_optimization_,
      min_branch_length_, max_iter_for_optimization_);
  branch_lengths_(op.gpcsp_) = branch_length;
}
