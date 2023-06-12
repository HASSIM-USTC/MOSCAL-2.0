#pragma once

#include "algebra.hpp"
#include "deom.hpp"
#include "linearalgebra.hpp"
#include "pulse.hpp"

#if defined(FERMI_QUAD)
#include "fermi_quad/fermi_quad.hpp"
#elif defined(FERMI_LINEAR)
#include "fermi/fermi.hpp"
#elif defined(BOSE_QUAD)
#include "bose_quad/bose_quad.hpp"
#elif defined(BOSE_LINEAR)
#include "bose/bose.hpp"
#elif defined(TEMPLATE)
#include "template/template.hpp"
#elif defined(BSM)
#include "bsm/bsm.hpp"
#elif defined(BSM_ACTION)
#include "bsm_action/bsm_action.hpp"
#elif defined(MD)
#include "md/md.hpp"
#endif

void A_assignments_B_p(int i) {}

template <typename Tc1, typename... Ts>
void A_assignments_B_p(int i, Tc1 &ddos, Tc1 &ddos1, Ts &...args) {
  ddos[i] = ddos1[i];
  A_assignments_B_p(i, args...);
}

template <typename... Ts>
void A_assignments_B(int i, int nsave6, Ts &...args) {
#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nsave6; i++) {
    A_assignments_B_p(i, args...);
  }
}

void clean_p(int i) {}

template <typename T, typename... Ts>
void clean_p(int i, T &arg, Ts &...args) {
  arg[i].setZero();
  clean_p(i, args...);
}

template <typename... Ts>
void clean(int i, int nsave6, Ts &...args) {
#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nsave6; i++) clean_p(i, args...);
}

template <class T>
void storage_eigen(T *syl, ComplexEigenSolver<MatrixXcd> eigensolver) {
  MatrixXcd V = eigensolver.eigenvectors();
  MatrixXcd VI = V.inverse();
  syl->S = eigensolver.eigenvalues();
  allocate_syl(syl->V, syl->VI, V, VI);
}

void construct_Mapping_filter(int nsave_1, int i, Trie *tree,
                              DEOM_DATA *d) {
#pragma omp for private(i) schedule(dynamic, 16)
  for (i = 0; i < nsave_1; i++) {
    construct_Mapping_filter_p(d, tree, i);
  }
}

void construct_Mapping_filter(nNNmat &ddos, int nsave_1, int i, Trie *tree,
                              DEOM_DATA *d) {
#pragma omp for private(i) schedule(dynamic, 16)
  for (i = 0; i < nsave_1; i++) {
    construct_Mapping_filter_p(ddos, d, tree, i);
  }
}

template <typename FUN>
void rem_cal_gen_filter(nNNmat &ddos, nNNmat &ddos1, int nsave_1, int i,
                        DEOM_DATA *d, FUN rem_cal_filter_p) {
#pragma omp for private(i) schedule(dynamic, 16)
  for (i = 0; i < nsave_1; i++) {
    rem_cal_filter_p(ddos, ddos1, d, i);
  }
}

void rem_cal(nNNmat &ddos, nNNmat &ddos1, int nsave_1, int i,
             DEOM_DATA *d) {
  rem_cal_gen_filter(ddos, ddos1, nsave_1, i, d, rem_cal_p);
}

void rem_cal_1(nNNmat &ddos, nNNmat &ddos1, int nsave_1, int i,
               DEOM_DATA *d) {
  rem_cal_gen_filter(ddos, ddos1, nsave_1, i, d, rem_cal_1_p);
}

void rem_cal_filter(nNNmat &ddos, nNNmat &ddos1, int nsave_1, int i,
                    DEOM_DATA *d) {
  rem_cal_gen_filter(ddos, ddos1, nsave_1, i, d, rem_cal_filter_p);
}

void rem_cal_1_filter(nNNmat &ddos, nNNmat &ddos1, int nsave_1, int i,
                      DEOM_DATA *d) {
  rem_cal_gen_filter(ddos, ddos1, nsave_1, i, d, rem_cal_1_filter_p);
}

void rem_cal_hei_filter(nNNmat &ddos, nNNmat &ddos1, int nsave_1, int i,
                        DEOM_DATA *d) {
  rem_cal_gen_filter(ddos, ddos1, nsave_1, i, d, rem_cal_hei_filter_p);
}

void rem_cal_1_hei_filter(nNNmat &ddos, nNNmat &ddos1, int nsave_1, int i,
                          DEOM_DATA *d) {
  rem_cal_gen_filter(ddos, ddos1, nsave_1, i, d, rem_cal_1_hei_filter_p);
}

void rem_cal_imag_filter(nNNmat &ddos, nNNmat &ddos1, int nsave_1, int i,
                         DEOM_DATA *d) {
  rem_cal_gen_filter(ddos, ddos1, nsave_1, i, d, rem_cal_imag_filter_p);
}

void rem_cal_cfw_filter(nNNmat &ddos, nNNmat &ddos1, int nsave_1, int i,
                        DEOM_DATA *d) {
  rem_cal_gen_filter(ddos, ddos1, nsave_1, i, d, rem_cal_cfw_filter_p);
}

template <typename FUN, typename... Ts>
void syl_cal_gen_filter(nNNmat &ddos, nNNmat &ddos1, int nsave_1, int i,
                        DEOM_DATA *d, DEOM_SYL *syl, FUN syl_cal_p,
                        Ts &...args) {
#pragma omp for private(i) schedule(dynamic, 16)
  for (i = 0; i < nsave_1; i++) {
    syl_cal_p(ddos, ddos1, args..., d, syl, i);
  }
}

void syl_cal_filter(nNNmat &ddos, nNNmat &ddos1, int nsave_1, int i,
                    DEOM_DATA *d, DEOM_SYL *syl) {
  syl_cal_gen_filter(ddos, ddos1, nsave_1, i, d, syl, syl_cal_filter_p);
}

void syl_cal_ddos_filter(nNNmat &ddos, nNNmat &ddos1, nNNmat &ddos2,
                         int nsave_1, int i, DEOM_DATA *d, DEOM_SYL *syl) {
  syl_cal_gen_filter(ddos, ddos1, nsave_1, i, d, syl,
                     syl_cal_ddos_filter_p, ddos2);
}

void syl_cal_1_filter(nNNmat &ddos, nNNmat &ddos1, nNNmat &ddos2,
                      int nsave_1, int i, DEOM_DATA *d, DEOM_SYL *syl) {
  syl_cal_gen_filter(ddos, ddos1, nsave_1, i, d, syl, syl_cal_1_filter_p,
                     ddos2);
}

void syl_cal_hei_filter(nNNmat &ddos, nNNmat &ddos1, nNNmat &ddos2,
                        int nsave_1, int i, DEOM_DATA *d, DEOM_SYL *syl) {
  syl_cal_gen_filter(ddos, ddos1, nsave_1, i, d, syl, syl_cal_hei_filter_p,
                     ddos2);
}

void syl_cal_1_hei_filter(nNNmat &ddos, nNNmat &ddos1, nNNmat &ddos2,
                          int nsave_1, int i, DEOM_DATA *d,
                          DEOM_SYL *syl) {
  syl_cal_gen_filter(ddos, ddos1, nsave_1, i, d, syl,
                     syl_cal_1_hei_filter_p, ddos2);
}

void syl_cal(nNNmat &ddos, nNNmat &ddos1, int i, DEOM_DATA *d,
             DEOM_SYL *syl) {
  syl_cal_gen_filter(ddos, ddos1, d->nddo, i, d, syl, syl_cal_p);
}

void syl_ddos_cal(nNNmat &ddos, nNNmat &ddos1, nNNmat &ddos2, int i,
                  DEOM_DATA *d, DEOM_SYL *syl) {
  syl_cal_gen_filter(ddos, ddos1, d->nddo, i, d, syl, syl_cal_ddos_p,
                     ddos2);
}

void syl_cal_1(nNNmat &ddos, nNNmat &ddos1, nNNmat &ddos2, int i,
               DEOM_DATA *d, DEOM_SYL *syl) {
  syl_cal_gen_filter(ddos, ddos1, d->nddo, i, d, syl, syl_cal_1_p, ddos2);
}

void rem_oprt(nNNmat &ddos, nNNmat &ddos1, int nsave_1, int i,
              DEOM_DATA *d, const char lcr) {
#pragma omp for private(i) schedule(dynamic, 16)
  for (i = 0; i < nsave_1; i++) {
    rem_oprt_p(ddos, ddos1, d, i, lcr);
  }
}

void rem_oprt_filter(nNNmat &ddos, nNNmat &ddos1, int nsave_1, int i,
                     DEOM_DATA *d, const char lcr) {
#pragma omp for private(i) schedule(dynamic, 16)
  for (i = 0; i < nsave_1; i++) {
    rem_oprt_filter_p(ddos, ddos1, d, i, lcr);
  }
}

template <typename FUN>
void filter(int nsave_1, int i, Trie *tree, DEOM_DATA *d, FUN filter_p) {
#pragma omp for private(i) schedule(dynamic, 16)
  for (i = 0; i < nsave_1; i++) {
    filter_p(d, d->ddos, d->keys, d->ddos1, d->keys1, tree, i);
  }
}

template <typename FUN>
void filter_2(int nsave_1, int i, Trie *tree, DEOM_DATA *d, FUN filter_p) {
#pragma omp for private(i) schedule(dynamic, 16)
  for (i = 0; i < nsave_1; i++)
    filter_p(d, d->ddos, d->ddos2, d->keys, d->ddos1, d->ddos3, d->keys1,
             tree, i);
}

void nsave_print(vector<int> &nsave, int nddo, int i) {
#pragma omp single
  {
    nsave[i] = nddo;
    printf("nddo:%12d\t", nsave[i]);
    fflush(stdout);
  }
}

void pulse_change_system(DEOM_DATA *d, const double t) {
  d->hamt = d->ham1 + pulse_et_ham1(t, d->pulse);
}