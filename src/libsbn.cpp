// Copyright 2019 Matsen group.
// libsbn is free software under the GPLv3; see LICENSE file for details.

#include "libsbn.hpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <random>
#include <string>

namespace py = pybind11;

class Matrix {
 public:
  Matrix(size_t rows, size_t cols) : m_rows(rows), m_cols(cols) {
    m_data = new double[rows * cols];
  }
  double *data() { return m_data; }
  size_t rows() const { return m_rows; }
  size_t cols() const { return m_cols; }

 private:
  size_t m_rows, m_cols;
  double *m_data;
};

class Vector {
 public:
  Vector(size_t size) : m_size(size) { m_data = new double[size]; }
  double *data() { return m_data; }
  size_t size() const { return m_size; }

 private:
  size_t m_size;
  double *m_data;
};

Vector make_vector() {
  std::random_device
      rd;  // Will be used to obtain a seed for the random number engine
  std::mt19937 gen(rd());  // Standard mersenne_twister_engine seeded with rd()
  std::uniform_real_distribution<> dis(0., 1.);
  Vector x(1000);
  for (size_t i = 0; i < x.size(); i++) {
    x.data()[i] = dis(gen);
  }
  return x;
}

double get00(Matrix m) { return m.data()[0]; }

PYBIND11_MODULE(sbn, m) {
  m.doc() = "libsbn bindings";
  py::class_<SBNInstance>(m, "instance")
      .def(py::init<const std::string &>())
      .def("tree_count", &SBNInstance::TreeCount)
      .def("read_newick_file", &SBNInstance::ReadNewickFile)
      .def("read_nexus_file", &SBNInstance::ReadNexusFile)
      .def("read_fasta_file", &SBNInstance::ReadFastaFile)
      .def("print_status", &SBNInstance::PrintStatus)
      .def("split_supports", &SBNInstance::SplitSupports)
      .def("make_beagle_instances", &SBNInstance::MakeBeagleInstances)
      .def("tree_log_likelihoods", &SBNInstance::TreeLogLikelihoods);
  m.def("f", &SBNInstance::f, "test");
  py::class_<Matrix>(m, "Matrix", py::buffer_protocol())
      .def_buffer([](Matrix &m) -> py::buffer_info {
        return py::buffer_info(
            m.data(),                                /* Pointer to buffer */
            sizeof(double),                          /* Size of one scalar */
            py::format_descriptor<double>::format(), /* Python struct-style
                                                       format descriptor */
            2,                                       /* Number of dimensions */
            {m.rows(), m.cols()},                    /* Buffer dimensions */
            {sizeof(double) * m.rows(), /* Strides (in bytes) for each index */
             sizeof(double)});
      });
  py::class_<Vector>(m, "Vector", py::buffer_protocol())
      .def_buffer([](Vector &v) -> py::buffer_info {
        return py::buffer_info(
            v.data(),                                /* Pointer to buffer */
            sizeof(double),                          /* Size of one scalar */
            py::format_descriptor<double>::format(), /* Python
                                                       struct-style format
                                                       descriptor */
            1,                                       /* Number of dimensions */
            {v.size()},                              /* Buffer dimensions */
            {sizeof(double)});
      });
  m.def("make_vector", &make_vector, "test");
  m.def("get00", &get00, "test");
}

