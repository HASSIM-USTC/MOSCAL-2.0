#pragma once
#include "deom.hpp"
#include "index.hpp"
#if defined(FERMI)
#include "aldebra-fermi.hpp"
#else // !FERMI
#include "aldebra-bose.hpp"
#endif // FERMI

bool is_valid(const MatrixNcd &d_ddo, const double ferr = 1e-10) {
  return d_ddo.squaredNorm() > ferr * ferr;
}

void print_is_hermite(const MatrixNcd &d_ddo, const double ferr = 1e-10) {
  printf((MatrixXcd(d_ddo) - MatrixXcd(d_ddo).adjoint()).squaredNorm() > ferr * ferr ? "NH" : "HE");
}

bool is_valid_max(const MatrixNcd &d_ddo, const double ferr = 1e-10) {
  for (int i = 0; i < NSYS; i++)
    for (int j = 0; j < NSYS; j++)
      if (abs(d_ddo.coeff(i, j)) > ferr)
        return true;
  return false;
}

bool is_valid(const MatrixNcd &d_ddo, const complex<double> absl, const double ferr = 1e-10) {
  return (abs(absl) == 0 ? 1 : abs(absl)) * d_ddo.squaredNorm() > ferr * ferr;
}

bool is_valid(const nNNmat &d_ddos, const int iado, const double ferr = 1e-10) {
  return d_ddos[iado].squaredNorm() > ferr * ferr;
}

complex<double> gen_absl(key_vec key, DEOM_DATA *d) {
  complex<double> absl = 2 * d->ham1_abs;
  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    if (key(i_nind)) {
      absl += abs(d->expn(i_nind));
      absl += 2 * d->qmd1_abs(i_nind) * d->coef_abs(i_nind);
    } else {
      absl += 2 * d->qmd1_abs(i_nind) * d->coef_abs(i_nind);
    }
  }
  return absl;
}

//        abs                 normal
//  max element of ddos       F-norm
// default is abs
bool is_valid(key_vec key, const MatrixNcd &d_ddo, DEOM_DATA *d) {
#if defined(NORMAL)
  complex<double> absl = gen_absl(key, d);
  return is_valid(d_ddo, absl, d->ferr);
#else  // ABS
  return is_valid_max(d_ddo, d->ferr);
#endif // NORMAL
}

int generate_hash_value_syl(const ArrayXi &key, const DEOM_DATA *d, const DEOM_SYL *syl) {
  int hash_value = 1;
  int sum = 0;
  for (int i = 0; i < syl->nind; i++) {
    sum += key(i);
    hash_value += d->comb_list(sum + i, i + 1);
  }
  return (int)hash_value;
}

void filter_p(DEOM_DATA *d, const nNNmat &ddos, const array_key_vec &keys, nNNmat &ddos1, array_key_vec &keys1, Trie *tree, int iado) {
  bool flag = is_valid(keys[iado], ddos[iado], d);
  if (flag) {
    int nddo_tmp;
#pragma omp critical
    {
      nddo_tmp = d->nddo;
      d->nddo++;
    }
    llint hash_val = gen_hash_value(keys[iado], d);
    tree->insert(hash_val, nddo_tmp);
    ddos1[nddo_tmp] = ddos[iado];
    keys1[nddo_tmp] = keys[iado];
  }
}

void filter_2_p(DEOM_DATA *d, const nNNmat &ddos, const nNNmat &ddos2, const array_key_vec &keys, nNNmat &ddos1, nNNmat &ddos3, array_key_vec &keys1, Trie *tree, int iado) {
  bool flag = is_valid(keys[iado], ddos[iado], d);
  if (flag) {
    int nddo_tmp;
#pragma omp critical
    {
      nddo_tmp = d->nddo;
      d->nddo++;
    }
    llint hash_val = gen_hash_value(keys[iado], d);
    tree->insert(hash_val, nddo_tmp);
    ddos1[nddo_tmp] = ddos[iado];
    ddos3[nddo_tmp] = ddos2[iado];
    keys1[nddo_tmp] = keys[iado];
  }
}

void filter_hei_p(DEOM_DATA *d, const nNNmat &ddos, const nNNmat &ddos2,
                  const array_key_vec &keys, nNNmat &ddos1, nNNmat &ddos3,
                  array_key_vec &keys1, Trie *tree, int iado) {
  bool flag = is_valid(keys[iado], ddos[iado], d);
  if (flag || ddos2[iado].squaredNorm() > d->ferr) {
    int nddo_tmp;
#pragma omp critical
    {
      nddo_tmp = d->nddo;
      d->nddo++;
    }
    llint hash_val = gen_hash_value(keys[iado], d);
    tree->insert(hash_val, nddo_tmp);
    ddos1[nddo_tmp] = ddos[iado];
    ddos3[nddo_tmp] = ddos2[iado];
    keys1[nddo_tmp] = keys[iado];
  }
}

void filter_zero_p(DEOM_DATA *d, const nNNmat &ddos, const array_key_vec &keys, nNNmat &ddos1, array_key_vec &keys1, Trie *tree, int iado) {
  if (is_valid(ddos[iado], 1e-15)) {
    int nddo_tmp;
#pragma omp critical
    {
      nddo_tmp = d->nddo;
      d->nddo++;
    }
    llint hash_val = gen_hash_value(keys[iado], d);
    tree->insert(hash_val, nddo_tmp);
    ddos1[nddo_tmp] = ddos[iado];
    keys1[nddo_tmp] = keys[iado];
  }
}

void set_ddos(const key_vec &key_swell, int iado, const int pos_get, const int hash_val, Trie *tree, DEOM_DATA *d) {
  int pos = tree->find(hash_val);
  if (pos == -1) {
    int nddo_tmp;
#pragma omp critical
    {
      nddo_tmp = d->nddo;
      d->nddo++;
    }
    const auto [rank_find, success] = tree->try_insert(hash_val, nddo_tmp);
    if (success) {
      d->g_index_list_pos[nddo_tmp](iado) = pos_get + 1;
      d->keys[nddo_tmp] = key_swell;
    } else {
      int nddo_tmp2 = tree->find(hash_val);
      d->g_index_list_pos[nddo_tmp2](iado) = pos_get + 1;
      d->keys[nddo_tmp].setZero();
    }
  } else {
    d->g_index_list_pos[pos](iado) = pos_get + 1;
  }
}
