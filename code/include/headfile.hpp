#pragma once
#define EIGEN_DONT_PARALLELIZE
// #define FMT_HEADER_ONLY

#if defined(FERMI_QUAD)
#define FERMI
#endif

#if defined(FERMI_LINEAR)
#define FERMI
#endif

#include <omp.h>
#include <openssl/md5.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <atomic>
#include <chrono>
#include <cmath>
#include <complex>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <thread>  // std::thread

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/SparseCore"
#include "index.hpp"
#include "nlohmann/json.hpp"

using Json = nlohmann::json;
using namespace std;
using namespace Eigen;

#if defined(STD)
#include "hashmap-std.hpp"
#else  // FOLLY
#include "hashmap-folly.hpp"
#endif  // STD

#if defined(SPARSE)
typedef SparseMatrix<complex<double>> MatrixNcd;
typedef SparseMatrix<complex<double>, Eigen::RowMajor> MatrixNcd_row;
#else
typedef Eigen::Matrix<complex<double>, NSYS, NSYS> MatrixNcd;
typedef Eigen::Matrix<complex<double>, NSYS, NSYS, Eigen::RowMajor>
    MatrixNcd_row;
#endif  // SPARSE

#if defined(FERMI)
typedef Eigen::Matrix<bool, 1, Eigen::Dynamic, Eigen::RowMajor> key_vec;
#else   // !FERMI
typedef VectorXi key_vec;
#endif  // FERMI

typedef Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> MatrixNb;
// max number of bytes used to store a pointer in 64 bit
typedef long long int llint;
typedef Eigen::Matrix<llint, Eigen::Dynamic, Eigen::Dynamic> llmat;
typedef std::vector<VectorXcd, Eigen::aligned_allocator<VectorXcd>> nXvec;
typedef std::vector<MatrixNcd, Eigen::aligned_allocator<MatrixNcd>> nNNmat;
typedef std::vector<MatrixNcd_row, Eigen::aligned_allocator<MatrixNcd_row>>
    nNNmat_row;
typedef std::vector<VectorXi, Eigen::aligned_allocator<VectorXi>>
    array_g_index;
typedef std::vector<key_vec, Eigen::aligned_allocator<key_vec>>
    array_key_vec;

void init_ddos_1_p(int nmax) {}

template <typename Tc1, typename... Ts>
void init_ddos_1_p(int nmax, Tc1 &ddos, Ts &...args) {
  ddos = nNNmat(nmax);
  init_ddos_1_p(nmax, args...);
}

void init_ddos_2_p(int i) {}

template <typename Tc1, typename... Ts>
void init_ddos_2_p(int i, Tc1 &ddos, Ts &...args) {
  ddos[i].resize(NSYS, NSYS);
  init_ddos_2_p(i, args...);
}

template <typename... Ts>
void init_ddos(int nmax, Ts &...args) {
  int i = 0;
  init_ddos_1_p(nmax, args...);

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nmax; i++) {
    init_ddos_2_p(i, args...);
  }
}

void init_nNNmat(int length) {}

template <typename Tc1, typename... Ts>
void init_nNNmat(int length, Tc1 &ddos, Ts &...args) {
  ddos = nNNmat(length);
  for (int i = 0; i < length; i++) ddos[i].resize(NSYS, NSYS);
  init_nNNmat(length, args...);
}

void init_nNNmat(MatrixNcd &mat) {
#ifdef SPARSE
  mat = mat.pruned(1);
  mat.makeCompressed();
#endif  // SPARSE matrix
}

void allocate_pulse(MatrixNcd &mat_out, MatrixXcd &mat_in) {
#ifdef SPARSE
  mat_out = (mat_in.sparseView(1));
  mat_out.makeCompressed();
#else   // dense matrix
  mat_out = mat_in;
#endif  // SPARSE matrix
}

void allocate_syl(MatrixNcd &sylU, MatrixNcd &sylUI, MatrixXcd &U,
                  MatrixXcd &UI) {
#ifdef SPARSE
  sylU = (U.sparseView(1));
  sylU.makeCompressed();
  sylUI = (UI.sparseView(1));
  sylUI.makeCompressed();
#else   // dense matrix
  sylU = U;
  sylUI = UI;
#endif  // SPARSE matrix
}

template <typename T>
void insert(MatrixNcd &mat, const int i, const int j, const T element) {
#ifdef SPARSE
  mat.insert(i, j) = element;
#else
  mat(i, j) = element;
#endif  // SPARSE
}