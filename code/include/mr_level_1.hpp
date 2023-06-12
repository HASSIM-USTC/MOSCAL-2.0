#pragma once

#include "deom.hpp"
#include "index.hpp"
#include "mr_level_2.hpp"
#include "print.hpp"
#include "rk.hpp"

#pragma omp declare reduction(+ : std::complex <double> : omp_out += omp_in)
complex<double> caculate_trace(DEOM_DATA *d, Trie *tree, Trie *tree_hei) {
  complex<double> trace = 0.0;
  int iado_;
#pragma omp parallel default(shared)
#pragma omp for private(iado_) reduction(+ : trace) schedule(dynamic, 16)
  for (iado_ = 0; iado_ < d->nddo; iado_++) {
    key_vec key(d->nind);
    key = d->keys[iado_];
    llint hash_val = gen_hash_value(key, d);
    int pos = tree->find(hash_val);
    int pos1 = tree_hei->find(hash_val);
    if ((pos1 != -1) && (pos == iado_)) {
      trace += MatrixXcd(d->ddos[pos] * d->ddos4[pos1]).trace();
    }
  }
  return trace;
}

complex<double> caculate_trace_filter(DEOM_DATA *d) {
  return caculate_trace(d, d->tree);
}

complex<double> caculate_trace_hei_filter(DEOM_DATA *d) {
  return caculate_trace(d, d->tree, d->tree_hei);
}

void init_nddo_clear_tree(DEOM_DATA *d) {
  d->nsave[6] = d->nddo;
  d->tree->clear();
  d->nddo = 0;
  d->nsave[0] = 0;
}

template <typename... Ts>
void check_if_converged(const CTRL *c, DEOM_DATA *d, Ts... args) {
  d->if_converged.trace = caculate_trace(d, args...);
  print(d->if_converged.trace, double(d->if_converged.flag_counter));
  if (abs((d->if_converged.trace_1 - d->if_converged.trace)) < c->ferr) {
    d->if_converged.flag_counter += 1;
    if (d->if_converged.flag_counter > c->step)
      d->if_converged.flag = false;
  } else {
    d->if_converged.flag_counter = 0;
  }
  d->if_converged.trace_1 = d->if_converged.trace;
}

void check_if_converged_filter(const CTRL *c, DEOM_DATA *d) {
  check_if_converged(c, d, d->tree);
}

void check_if_converged_hei_filter(const CTRL *c, DEOM_DATA *d) {
  check_if_converged(c, d, d->tree, d->tree_hei);
}

void renew_check_flag(DEOM_DATA *d) {
  d->if_converged.flag_counter = 0;
  d->if_converged.flag = true;
}

void generate_syl_omega(DEOM_DATA *d, DEOM_SYL *syl, double w) {
  // This is S and B in $AX - XB = C$ in syl_cal_filter
  // But is R and A in $AX - XB = C$ in syl_cal_Hei_filter
  ComplexEigenSolver<MatrixXcd> eigensolver(
      (MatrixXcd)(d->i * (d->hamt + w * d->I) - syl->OMG * d->I));
  storage_eigen(syl, eigensolver);
}

void Init_ddos_hei(DEOM_DATA *d, nNNmat &ddos) {
  for (int k = 0; k < d->nmod; k++) ddos[0] += d->dipole1.sdip[k];
  d->tree->try_insert(1, 0);

  key_vec key_swell(d->nind);
  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    llint hash_val = generate_key_plus(d->zerokey, key_swell, i_nind, d);
    int nddo_tmp = d->nddo;
    d->nddo++;
    const auto [rank_find, success] =
        d->tree->try_insert(hash_val, nddo_tmp);
    if (success) {
      d->keys[rank_find] = key_swell;
      ddos[rank_find] = d->dipole1.qmdta_l[i_nind];
    }
  }
}

template <typename FUN>
void oprt(const char lcr, int i, DEOM_DATA *d, FUN rem_oprt) {
#pragma omp parallel default(shared)
  {
    rem_oprt(d->ddos1, d->ddos, d->nddo, i, d, lcr);
    A_assignments_B(i, d->nddo, d->ddos, d->ddos1);
  }
}

template <typename FUN>
void oprt(const char lcr, int i, Trie *tree, DEOM_DATA *d,
          FUN rem_oprt_filter) {
  init_nddo_clear_tree(d);

#pragma omp parallel default(shared)
  {
    clean(i, d->nsave[6], d->g_index_list_pos, d->ddos1);
    filter(d->nsave[6], i, tree, d, filter_p);
    A_assignments_B(i, d->nsave[6], d->ddos, d->ddos1, d->keys, d->keys1);
    nsave_print(d->nsave, d->nddo, 1);

    construct_Mapping_filter(d->nsave[1], i, tree, d);
    nsave_print(d->nsave, d->nddo, 2);
    rem_oprt_filter(d->ddos1, d->ddos, d->nsave[2], i, d, lcr);
    A_assignments_B(i, d->nsave[2], d->ddos, d->ddos1, d->ddos2, d->ddos1);
  }
}

template <typename FUN>
void oprt_filter(const char lcr, int i, DEOM_DATA *d, FUN rem_oprt) {
  oprt(lcr, i, d->tree, d, rem_oprt);
}

void allocator_system(int i, DEOM_DATA *d, const char lcr) {
  init_ddos(d->nmax, d->ddos4);
  oprt(lcr, i, d->tree_hei, d, rem_oprt_filter);

#pragma omp parallel default(shared)
  {
    A_assignments_B(i, d->nsave[2], d->ddos4, d->ddos1);
    clean(i, d->nmax, d->ddos, d->ddos1, d->ddos2, d->ddos3, d->keys,
          d->keys1, d->g_index_list_pos);
  }

  d->nddo = 1;
}

void allocator_bath(int i, DEOM_DATA *d, const char lcr) {
  init_ddos(d->nmax, d->ddos4);
  oprt(lcr, i, d->tree_hei, d, rem_oprt_filter);

#pragma omp parallel default(shared)
  {
    A_assignments_B(i, d->nsave[2], d->ddos4, d->ddos1);
    clean(i, d->nmax, d->ddos, d->ddos1, d->ddos2, d->ddos3, d->keys,
          d->keys1, d->g_index_list_pos);
  }

  d->nddo = 1;
}

void gen_keys(DEOM_DATA *d) {
  int i, filter = d->lmax + 1, nddo_bef = 0;
  d->nddo = 1;

#pragma omp parallel default(shared)
  for (int ii = 0; ii < filter; ii++) {
#pragma omp for private(i) schedule(dynamic, 16)
    for (i = nddo_bef; i < d->nddo; i++) construct_Mapping_p(d, i);

#pragma omp single
    {
      nddo_bef = d->nddo;
#if defined(FERMI)
      d->nddo = 0;
      for (i = 0; i <= ii; i++)
        if (i <= d->lmax) {
          d->nddo += d->comb_list(d->nind, i);
        }
#else
      if (ii <= d->lmax) {
        d->nddo = d->comb_list(d->nind + ii, d->nind);
      } else {
        d->nddo = d->comb_list(d->nind + d->lmax, d->nind);
      }
#endif  // FERMI
      printf("nddo:%d\n", d->nddo);
    }
  }
}

void gen_syl_p(DEOM_DATA *d, DEOM_SYL *syl, int iado) {
  ArrayXi key_syl(syl->nind);
  key_syl = syl->keys.row(iado);
  complex<double> gams = 0;
  for (int i_nind = 0; i_nind < syl->nind; i_nind++)
    gams += double(key_syl(i_nind)) * syl->expn(i_nind);

  MatrixXcd hlft = d->i * d->hamt + gams * d->I;
  ComplexEigenSolver<MatrixXcd> eigensolverA(hlft);
  auto R = eigensolverA.eigenvalues();
  MatrixXcd U = eigensolverA.eigenvectors();
  MatrixXcd UI = U.inverse();
  syl->R[iado] = R;
  allocate_syl(syl->U[iado], syl->UI[iado], U, UI);
}

void construct_Mapping_syl_p(DEOM_DATA *d, DEOM_SYL *syl, int iado) {
  ArrayXi key(syl->nind), key_swell(syl->nind);
  key = syl->keys.row(iado);
  for (int i_nind = 0; i_nind < syl->nind; i_nind++) {
    if (key.sum() < d->lmax) {
      key_swell = key;
      key_swell(i_nind) += 1;
      llint hash_val = generate_hash_value_syl(key_swell, d, syl) - 1;
      if (hash_val != -1) {
        syl->keys.row(hash_val) = key_swell;
      }
    }
  }
}

void gen_syl(DEOM_DATA *d, DEOM_SYL *syl) {
  int i, filter = d->lmax + 1, nddo_bef = 0;
  syl->nddo = 1;
#pragma omp parallel default(shared)
  {
    for (int ii = 0; ii < filter; ii++) {
#pragma omp for private(i) schedule(dynamic, 16)
      for (i = nddo_bef; i < syl->nddo; i++) {
        construct_Mapping_syl_p(d, syl, i);
      }

#pragma omp single
      {
        nddo_bef = syl->nddo;
        syl->nddo = d->comb_list(syl->nind + ii, ii);
        printf("nddo:%d\tnind:%d\tlmax:%d\n", syl->nddo, syl->nind,
               d->lmax);
      }
    }

#pragma omp for private(i) schedule(dynamic, 16)
    for (i = 0; i < syl->nddo; i++) {
      gen_syl_p(d, syl, i);
    }
  }
}

template <typename FUN>
void sc2_ddos_cal(DEOM_DATA *d, DEOM_SYL *syl, int i, FUN syl_cal_1) {
#pragma omp parallel default(shared)
  {
    syl_cal_1(d->ddos1, d->ddos, d->ddos2, i, d, syl);
    A_assignments_B(i, d->nddo, d->ddos, d->ddos1);
  }
}

template <typename FUN>
void sc2_cal(DEOM_DATA *d, DEOM_SYL *syl, int i, FUN syl_cal_1) {
#pragma omp parallel default(shared)
  {
    syl_cal_1(d->ddos1, d->ddos, i, d, syl);
    A_assignments_B(i, d->nddo, d->ddos, d->ddos1);
  }
}

template <typename FUN>
void sci_calculator_equal_ddos_filter_1(DEOM_DATA *d, int i,
                                        FUN filter_2_p) {
#pragma omp parallel default(shared)
  {
    clean(i, d->nsave[6], d->g_index_list_pos, d->ddos1, d->ddos3);
    filter_2(d->nsave[6], i, d->tree, d, filter_2_p);
    A_assignments_B(i, d->nsave[6], d->ddos, d->ddos1, d->ddos2, d->ddos3,
                    d->keys, d->keys1);
    nsave_print(d->nsave, d->nddo, 1);
  }
}

template <typename FUN>
void sc2_cal_filter_1(DEOM_DATA *d, int i, FUN filter_p) {
#pragma omp parallel default(shared)
  {
    clean(i, d->nsave[6], d->g_index_list_pos, d->ddos1);
    filter(d->nsave[6], i, d->tree, d, filter_p);
    A_assignments_B(i, d->nsave[6], d->ddos, d->ddos1, d->keys, d->keys1);
    nsave_print(d->nsave, d->nddo, 1);
  }
}

template <typename FUN>
void sci_calculator_equal_ddos_filter_2(DEOM_DATA *d, DEOM_SYL *syl, int i,
                                        FUN syl_cal) {
#pragma omp parallel default(shared)
  {
    construct_Mapping_filter(d->nsave[1], i, d->tree, d);
    nsave_print(d->nsave, d->nddo, 2);
    syl_cal(d->ddos1, d->ddos, d->ddos2, d->nsave[2], i, d, syl);
    A_assignments_B(i, d->nsave[2], d->ddos, d->ddos1);
  }
}

template <typename FUN>
void sc2_cal_filter_2(DEOM_DATA *d, DEOM_SYL *syl, int i, FUN syl_cal) {
#pragma omp parallel default(shared)
  {
    construct_Mapping_filter(d->nsave[1], i, d->tree, d);
    nsave_print(d->nsave, d->nddo, 2);
    syl_cal_filter(d->ddos1, d->ddos, d->nsave[2], i, d, syl);
    A_assignments_B(i, d->nsave[2], d->ddos, d->ddos1);
  }
}