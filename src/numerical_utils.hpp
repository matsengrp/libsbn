// Copyright 2019 libsbn project contributors.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#ifndef SRC_NUMERICAL_UTILS_HPP_
#define SRC_NUMERICAL_UTILS_HPP_

#include <limits.h>

#include "eigen_sugar.hpp"
#include "sugar.hpp"

constexpr double DOUBLE_INF = std::numeric_limits<double>::infinity();
constexpr double DOUBLE_NEG_INF = -std::numeric_limits<double>::infinity();
constexpr double EPS = std::numeric_limits<double>::epsilon();
constexpr double LOG_EPS = log(EPS);

namespace NumericalUtils {

// Return log(exp(x) + exp(y)).
constexpr double LogAdd(double x, double y) {
  // See:
  // https://github.com/alexandrebouchard/bayonet/blob/master/src/main/java/bayonet/math/NumericalUtils.java#L59

  // Make x the max.
  if (y > x) {
    double temp = x;
    x = y;
    y = temp;
  }
  if (x == DOUBLE_NEG_INF) {
    return x;
  }
  double negDiff = y - x;
  if (negDiff < LOG_EPS) {
    return x;
  }
  return x + log(1.0 + exp(negDiff));
}

// Return log(sum_i exp(vec(i))).
double LogSum(const EigenVectorXdRef vec);
// Normalize the entries of vec: vec(i) = vec(i) - LogSum(vec).
void ProbabilityNormalizeInLog(EigenVectorXdRef vec);
// Exponentiate vec in place: vec(i) = exp(vec(i))
void Exponentiate(EigenVectorXdRef vec);

}  // namespace NumericalUtils

#ifdef DOCTEST_LIBRARY_INCLUDED
TEST_CASE("NumericalUtils") {
  double log_x = log(2);
  double log_y = log(3);
  double log_sum = NumericalUtils::LogAdd(log_x, log_y);
  CHECK_LT(fabs(log_sum - 1.609438), 1e-5);

  EigenVectorXd log_vec(10);
  double log_sum2 = DOUBLE_NEG_INF;
  for (size_t i = 0; i < log_vec.size(); i++) {
    log_vec(i) = log(i + 1);
    log_sum2 = NumericalUtils::LogAdd(log_sum2, log_vec(i));
  }
  log_sum = NumericalUtils::LogSum(log_vec);
  CHECK_LT(fabs(log_sum - 4.007333), 1e-5);
  CHECK_LT(fabs(log_sum2 - 4.007333), 1e-5);

  NumericalUtils::ProbabilityNormalizeInLog(log_vec);
  for (size_t i = 0; i < log_vec.size(); i++) {
    CHECK_LT(fabs(log_vec(i) - (log(i + 1) - log_sum)), 1e-5);
  }

  NumericalUtils::Exponentiate(log_vec);
  double sum = 0.0;
  for (size_t i = 0; i < log_vec.size(); i++) {
    sum += log_vec(i);
  }
  CHECK_LT(fabs(sum - 1), 1e-5);
}
#endif  // DOCTEST_LIBRARY_INCLUDED

#endif /* SRC_NUMERICAL_UTILS_HPP_ */
