// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.
//
// Perform training of an SBN based on a sample of trees.
//
// We assume that readers are familiar with how the sbn_parameters_ vector is laid out:
// first probabilities of rootsplits, then conditional probabilities of PCSSs.

#include "sbn_probability.hpp"
#include <numeric>
#include "sbn_maps.hpp"

// Increment all entries from an index vector by a value.
void IncrementBy(EigenVectorXdRef vec, const SizeVector& indices, double value) {
  for (const auto& idx : indices) {
    vec[idx] += value;
  }
}

// Increment all entries from an index vector vector by a value.
void IncrementBy(EigenVectorXdRef vec, const SizeVectorVector& index_vector_vector,
                 double value) {
  for (const auto& indices : index_vector_vector) {
    IncrementBy(vec, indices, value);
  }
}

// Increment all entries from an index vector of length k by a value drawn from a vector
// of values of length k.
void IncrementBy(EigenVectorXdRef vec, const SizeVector& indices,
                 const EigenConstVectorXdRef values) {
  Assert(indices.size() == values.size(),
         "Indices and values don't have matching size.");
  for (size_t i = 0; i < values.size(); ++i) {
    vec[indices[i]] += values[i];
  }
}

// Repeat the previous increment operation across a vector of index vectors with a fixed
// vector of values.
void IncrementBy(EigenVectorXdRef vec, const SizeVectorVector& index_vector_vector,
                 const EigenConstVectorXdRef values) {
  Assert(index_vector_vector.size() == values.size(),
         "Indices and values don't have matching size.");
  for (size_t i = 0; i < values.size(); ++i) {
    IncrementBy(vec, index_vector_vector[i], values[i]);
  }
}

// I know I'm supposed to prefer
// return std::accumulate(indices.cbegin(), indices.cend(), starting_value,
//                        [&vec](const double& subproduct, const size_t& idx) {
//                          return subproduct * vec[idx];
//                        });
// but the version below is so clear that it doesn't need any comments...
double ProductOf(const EigenConstVectorXdRef vec, const SizeVector& indices,
                 const double starting_value) {
  double result = starting_value;
  for (const auto& idx : indices) {
    result *= vec[idx];
  }
  return result;
}

void ProbabilityNormalizeRange(EigenVectorXdRef vec, std::pair<size_t, size_t> range) {
  auto [start_idx, end_idx] = range;
  auto segment = vec.segment(start_idx, end_idx - start_idx);
  segment /= segment.sum();
}

// We assume that vec is laid out like sbn_parameters (see top).
void ProbabilityNormalizeParams(EigenVectorXdRef vec, size_t rootsplit_count,
                                const BitsetSizePairMap& parent_to_range) {
  ProbabilityNormalizeRange(vec, {0, rootsplit_count});
  for (const auto& [_, range] : parent_to_range) {
    ProbabilityNormalizeRange(vec, range);
  }
}

IndexerRepresentationCounter SBNProbability::IndexerRepresentationCounterOf(
    const BitsetSizeMap& indexer, const Node::TopologyCounter& topology_counter) {
  IndexerRepresentationCounter counter;
  counter.reserve(topology_counter.size());
  for (const auto& [topology, topology_count] : topology_counter) {
    counter.push_back(
        {SBNMaps::IndexerRepresentationOf(indexer, topology), topology_count});
  }
  return counter;
}

// Set the provided counts vector to be the counts of the rootsplits and PCSSs provided
// in the input.
void AccumulateCounts(
    EigenVectorXdRef counts,
    const IndexerRepresentationCounter& indexer_representation_counter,
    size_t rootsplit_count, const BitsetSizePairMap& parent_to_range) {
  counts.setZero();
  for (const auto& [indexer_representation, int_topology_count] :
       indexer_representation_counter) {
    const auto& [rootsplits, pcss_vector_vector] = indexer_representation;
    const auto topology_count = static_cast<double>(int_topology_count);
    IncrementBy(counts, rootsplits, topology_count);
    IncrementBy(counts, pcss_vector_vector, topology_count);
  }
}

void SBNProbability::SimpleAverage(
    EigenVectorXdRef sbn_parameters,
    const IndexerRepresentationCounter& indexer_representation_counter,
    size_t rootsplit_count, const BitsetSizePairMap& parent_to_range) {
  AccumulateCounts(sbn_parameters, indexer_representation_counter, rootsplit_count,
                   parent_to_range);
  ProbabilityNormalizeParams(sbn_parameters, rootsplit_count, parent_to_range);
}

void SBNProbability::ExpectationMaximization(
    EigenVectorXdRef sbn_parameters,
    const IndexerRepresentationCounter& indexer_representation_counter,
    size_t rootsplit_count, const BitsetSizePairMap& parent_to_range, double alpha,
    size_t em_loop_count) {
  // This vector holds the \bar{m} vectors (described in the 2018 NeurIPS paper).
  // They are packed into a single vector as sbn_parameters is.
  EigenVectorXd m_bar(sbn_parameters.size());
  // The q weight of a rootsplit is the probability of each rooting given the current
  // SBN parameters.
  Assert(!indexer_representation_counter.empty(),
         "Empty indexer_representation_counter.");
  auto edge_count = indexer_representation_counter[0].first.first.size();
  EigenVectorXd q_weights(edge_count);
  // This vector holds the \tilde{m} vectors (described in the 2018 NeurIPS paper),
  // which is the counts vector before normalization to get the SimpleAverage estimate.
  EigenVectorXd m_tilde(sbn_parameters.size());
  AccumulateCounts(m_tilde, indexer_representation_counter, rootsplit_count,
                   parent_to_range);
  // The normalized version of m_tilde is the SA estimate, which is our starting point.
  sbn_parameters = m_tilde;
  ProbabilityNormalizeParams(sbn_parameters, rootsplit_count, parent_to_range);
  // Do the specified number of EM loops.
  for (size_t em_idx = 0; em_idx < em_loop_count; ++em_idx) {
    m_bar.setZero();
    // Loop over topologies (as manifested by their indexer representations).
    for (const auto& [indexer_representation, int_topology_count] :
         indexer_representation_counter) {
      const auto& [rootsplits, pcss_vector_vector] = indexer_representation;
      const auto topology_count = static_cast<double>(int_topology_count);
      // Calculate the q weights for this topology.
      q_weights.setZero();
      Assert(rootsplits.size() == edge_count,
             "Rootsplit length not equal to edge_count.");
      Assert(pcss_vector_vector.size() == edge_count,
             "PCSSs length not equal to edge_count.");
      // Loop over the various rooting positions of this topology.
      for (size_t rooting_position = 0; rooting_position < edge_count;
           ++rooting_position) {
        const size_t& rootsplit = rootsplits[rooting_position];
        const SizeVector& pcss_vector = pcss_vector_vector[rooting_position];
        // Calculate the SBN probability of this topology rooted at this position.
        q_weights[rooting_position] =
            ProductOf(sbn_parameters, pcss_vector, sbn_parameters[rootsplit]);
      }
      q_weights /= q_weights.sum();
      // Now increment the new SBN parameters by the q-weighted counts.
      q_weights *= topology_count;
      IncrementBy(m_bar, rootsplits, q_weights);
      IncrementBy(m_bar, pcss_vector_vector, q_weights);
    }
    sbn_parameters = m_bar + alpha * m_tilde;
    ProbabilityNormalizeParams(sbn_parameters, rootsplit_count, parent_to_range);
  }
}

double SBNProbability::ProbabilityOf(
    const EigenConstVectorXdRef sbn_parameters,
    const IndexerRepresentation& indexer_representation) {
  const auto& [rootsplits, pcss_vector_vector] = indexer_representation;
  auto single_rooting_probability = [&sbn_parameters](const size_t rootsplit,
                                                      const SizeVector pcss_vector) {
    return ProductOf(sbn_parameters, pcss_vector, sbn_parameters[rootsplit]);
  };
  return std::inner_product(
      rootsplits.cbegin(), rootsplits.cend(),  // First vector.
      pcss_vector_vector.cbegin(),             // Second vector.
      0.,                                      // Starting value.
      std::plus<>(),                           // "Reduce" using +.
      single_rooting_probability);             // How to combine pairs of elements.
}

EigenVectorXd SBNProbability::ProbabilityOf(
    const EigenConstVectorXdRef sbn_parameters,
    const std::vector<IndexerRepresentation>& indexer_representations) {
  const size_t tree_count = indexer_representations.size();
  EigenVectorXd results(tree_count);
  for (size_t tree_idx = 0; tree_idx < tree_count; ++tree_idx) {
    results[tree_idx] =
        ProbabilityOf(sbn_parameters, indexer_representations[tree_idx]);
  }
  return results;
}
