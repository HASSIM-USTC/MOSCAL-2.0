#pragma once

#include <Eigen/Eigenvalues>

#include "algebra.hpp"
#include "deom.hpp"

void gen_syl_key(ArrayXi &key_syl, key_vec &key, const DEOM_SYL *syl) {
  for (int i_nind = 0; i_nind < syl->nind; i_nind++) {
    key_syl(i_nind) = 0;
    for (int i = 0; i < syl->lwsg; i++)
      if (key(syl->twsg(i_nind, i))) key_syl(i_nind) += 1;
  }
}

#if defined(SPARSE)
template <typename T1, typename T2>
void syl_gen_Y(MatrixNcd &D, T1 &R, T2 &S, MatrixNcd &Y) {
  for (int k = 0; k < D.outerSize(); ++k) {
    for (MatrixNcd::InnerIterator it(D, k); it; ++it) {
      int i = it.row();
      int j = it.col();
      Y.insert(i, j) = it.value() / (R(i) - S(j));
    }
  }
}

void syl_complex(MatrixNcd &X, MatrixNcd &A, MatrixNcd &B, MatrixNcd &C) {
  ComplexEigenSolver<MatrixXcd> eigensolverA(A);
  ComplexEigenSolver<MatrixXcd> eigensolverB(B);
  auto R = eigensolverA.eigenvalues();
  auto S = eigensolverB.eigenvalues();
  MatrixNcd U = (eigensolverA.eigenvectors()).sparseView(1);
  MatrixNcd V = (eigensolverB.eigenvectors()).sparseView(1);
  MatrixNcd UI = (eigensolverA.eigenvectors().inverse()).sparseView(1);
  MatrixNcd VI = (eigensolverB.eigenvectors().inverse()).sparseView(1);
  MatrixNcd Y = MatrixNcd(NSYS, NSYS);
  Y.setZero();
  MatrixNcd D = ((UI * C).pruned(1) * V).pruned(1);
  syl_gen_Y(D, R, S, Y);
  X = (U * Y * VI).pruned(1);
}

#else   // !SPARSE

template <typename T1, typename T2>
void syl_gen_Y(MatrixNcd &D, T1 &R, T2 &S, MatrixNcd &Y) {
  for (int i = 0; i < NSYS; i++)
    for (int j = 0; j < NSYS; j++) Y(i, j) = D(i, j) / (R(i) - S(j));
}

void syl_complex(MatrixNcd &X, MatrixNcd &A, MatrixNcd &B, MatrixNcd &C) {
  ComplexEigenSolver<MatrixNcd> eigensolverA(A);
  ComplexEigenSolver<MatrixNcd> eigensolverB(B);
  auto R = eigensolverA.eigenvalues();
  auto S = eigensolverB.eigenvalues();
  MatrixNcd U = eigensolverA.eigenvectors();
  MatrixNcd V = eigensolverB.eigenvectors();
  MatrixNcd Y = MatrixNcd::Zero(NSYS, NSYS);
  MatrixNcd D = U.inverse() * C * V;
  syl_gen_Y(D, R, S, Y);
  X = (U * Y * V.inverse());
}
#endif  // SPARSE

void syl_complex_storge(MatrixNcd &X, const DEOM_SYL *syl, int pos0,
                        double ferr, MatrixNcd &C) {
  MatrixNcd Y = MatrixNcd(NSYS, NSYS);
  Y.setZero();
  MatrixNcd D = syl->UI[pos0] * C * syl->V;
  syl_gen_Y(D, syl->R[pos0], syl->S, Y);
  X = syl->U[pos0] * Y * syl->VI;
  init_nNNmat(X);
}

void syl_complex_Hei_storge(MatrixNcd &X, const DEOM_SYL *syl, int pos0,
                            double ferr, MatrixNcd &C) {
  MatrixNcd Y = MatrixNcd(NSYS, NSYS);
  Y.setZero();
  MatrixNcd D = syl->VI * C * syl->U[pos0];
  syl_gen_Y(D, syl->S, syl->R[pos0], Y);
  X = syl->V * Y * syl->UI[pos0];
  init_nNNmat(X);
}
