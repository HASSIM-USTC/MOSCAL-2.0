#pragma once

#include "algebra.hpp"
#include "deom.hpp"

void construct_Mapping_filter_p_no_add(DEOM_DATA *d, Trie *tree, int iado) {
  key_vec key(d->keys[iado]), key_swell(d->nind);
  llint hash_val = gen_hash_value(key, d);
  int pos_get = tree->find(hash_val);
  d->g_index_list_pos[pos_get](0) = pos_get + 1;

  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    if (key(i_nind)) {
      hash_val = generate_key_minus(key, key_swell, i_nind, d);
      int pos = tree->find(hash_val);
      if (pos != -1)
        d->g_index_list_pos[pos](i_nind + 1) = pos_get + 1;
    } else {
      if (tier(key, d) < d->lmax) {
        hash_val = generate_key_plus(key, key_swell, i_nind, d);
        int pos = tree->find(hash_val);
        if (pos != -1)
          d->g_index_list_pos[pos](i_nind + 1) = pos_get + 1;
      }
    }
  }
}

void construct_Mapping_filter_p(nNNmat &ddos, DEOM_DATA *d, Trie *tree, int iado) {
  key_vec key(d->nind);
  key_vec key_swell(d->nind);
  llint hash_val;
  MatrixNcd addo;
  addo.resize(NSYS, NSYS);
  key = d->keys[iado];
  hash_val = generate_hash_value(key, d);

  complex<double> absl = d->ham1_abs;
  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    if (key(i_nind)) {
      absl += abs(d->expn(i_nind));
      absl += d->qmd1_abs(i_nind);
    } else {
      absl += d->qmd1_abs(i_nind) * d->coef_abs(i_nind);
    }
  }

  if (is_valid(addo, absl, d->ferr) || hash_val == 1) {
    int pos_get = tree->find(hash_val);
    d->g_index_list_pos[pos_get](0) = pos_get + 1;
    for (int i_nind = 0; i_nind < d->nind; i_nind++) {
      if (key(i_nind)) {
        hash_val = generate_key_minus(key, key_swell, i_nind, d);
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
            d->g_index_list_pos[nddo_tmp](i_nind + 1) = pos_get + 1;
            d->keys[nddo_tmp] = key_swell;
          } else {
            pos = tree->find(hash_val);
            d->g_index_list_pos[pos](i_nind + 1) = pos_get + 1;
            d->keys[nddo_tmp].setZero();
          }
        } else {
          d->g_index_list_pos[pos](i_nind + 1) = pos_get + 1;
        }
      } else {
        if (tier(key, d) < d->lmax) {
          hash_val = generate_key_plus(key, key_swell, i_nind, d);
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
              d->g_index_list_pos[nddo_tmp](i_nind + 1) = pos_get + 1;
              d->keys[nddo_tmp] = key_swell;
            } else {
              pos = tree->find(hash_val);
              d->g_index_list_pos[pos](i_nind + 1) = pos_get + 1;
              d->keys[nddo_tmp].setZero();
            }
          } else {
            d->g_index_list_pos[pos](i_nind + 1) = pos_get + 1;
          }
        }
      }
    }
  }
}

void filter_p(DEOM_DATA *d, const nNNmat &ddos, const array_key_vec &keys, nNNmat &ddos1, array_key_vec &keys1, Trie *tree, int iado) {
  MatrixNcd addo, addo_fil;
  addo.resize(NSYS, NSYS);
  addo_fil.resize(NSYS, NSYS);
  key_vec key(d->nind);
  key = keys[iado];
  addo = ddos[iado];
  key_vec parity_key(d->nind);
  complex<double> gamma = 0.0;

  bool parity = true;
  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    if (key(i_nind)) {
      parity ^= true;
      gamma += d->expn(i_nind);
    }
    parity_key(i_nind) = parity;
  }
  double parity_p_n = parity ? -1.0 : 1.0;

  addo_fil = d->i * (d->ham1 * addo - addo * d->ham1) + gamma * addo;
  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    if (key(i_nind)) {
      addo_fil += (parity_key(i_nind) ? 1.0 : -1.0) * d->i * d->coef_abs(i_nind) * (d->qmdta[i_nind] * addo - addo * d->qmdta[i_nind]);
    } else {
      int mod = d->nmodmaxLabel(i_nind);
      for (int ii = 0; ii < d->nmodmax; ii++) {
        addo_fil += (parity_key(i_nind) ? -1.0 : 1.0) * d->i / d->coef_abs(i_nind) * ((double)parity_p_n * d->coef_lft[i_nind](mod, ii) * d->qmdtc[i_nind] * addo + d->coef_rht[i_nind](mod, ii) * addo * d->qmdtc[i_nind]);
      }
    }
  }
  if (addo_fil.squaredNorm() > d->ferr) {
    int nddo_tmp;
#pragma omp critical
    {
      nddo_tmp = d->nddo;
      d->nddo++;
    }
    llint hash_val = generate_hash_value(key, d);
    tree->insert(hash_val, nddo_tmp);
    ddos1[nddo_tmp] = addo;
    keys1[nddo_tmp] = key;
  }
}

void filter_2_p(DEOM_DATA *d, const nNNmat &ddos, const nNNmat &ddos2, const array_key_vec &keys, nNNmat &ddos1, nNNmat &ddos3, array_key_vec &keys1, Trie *tree, int iado) {
  key_vec key(d->nind);
  MatrixNcd addo, addo1, addo_fil;
  addo.resize(NSYS, NSYS);
  addo1.resize(NSYS, NSYS);
  addo_fil.resize(NSYS, NSYS);
  key = keys[iado];
  addo = ddos[iado];
  addo1 = ddos2[iado];
  key_vec parity_key(d->nind);
  complex<double> gamma = 0.0;

  bool parity = true;
  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    if (key(i_nind)) {
      parity ^= true;
      gamma += d->expn(i_nind);
    }
    parity_key(i_nind) = parity;
  }
  double parity_p_n = parity ? -1.0 : 1.0;

  addo_fil = d->i * (d->ham1 * addo - addo * d->ham1) + (gamma * addo);
  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    if (key(i_nind)) {
      addo_fil += (parity_key(i_nind) ? 1.0 : -1.0) * d->i * d->coef_abs(i_nind) * (d->qmdta[i_nind] * addo - addo * d->qmdta[i_nind]);
    } else {
      int mod = d->nmodmaxLabel(i_nind);
      for (int ii = 0; ii < d->nmodmax; ii++) {
        addo_fil += (parity_key(i_nind) ? -1.0 : 1.0) * d->i / d->coef_abs(i_nind) * ((double)parity_p_n * d->coef_lft[i_nind](mod, ii) * d->qmdtc[i_nind] * addo + d->coef_rht[i_nind](mod, ii) * addo * d->qmdtc[i_nind]);
      }
    }
  }

  if (addo_fil.squaredNorm() > d->ferr) {
    int nddo_tmp;
#pragma omp critical
    {
      nddo_tmp = d->nddo;
      d->nddo++;
    }
    llint hash_val = generate_hash_value(key, d);
    tree->insert(hash_val, nddo_tmp);
    ddos1[nddo_tmp] = addo;
    ddos3[nddo_tmp] = addo1;
    keys1[nddo_tmp] = key;
  }
}

void filter_hei_p(DEOM_DATA *d, const nNNmat &ddos, const nNNmat &ddos2, const array_key_vec &keys, nNNmat &ddos1, nNNmat &ddos3, array_key_vec &keys1, Trie *tree, int iado) {
  key_vec key(d->nind);
  MatrixNcd addo, addo1, addo_fil;
  addo.resize(NSYS, NSYS);
  addo1.resize(NSYS, NSYS);
  addo_fil.resize(NSYS, NSYS);
  key = keys[iado];
  addo = ddos[iado];
  addo1 = ddos2[iado];
  key_vec parity_key(d->nind);
  complex<double> gamma = 0.0;

  bool parity = true;
  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    if (key(i_nind)) {
      parity ^= true;
      gamma += d->expn(i_nind);
    }
    parity_key(i_nind) = parity;
  }
  double parity_p_n = parity ? -1.0 : 1.0;
  addo_fil = d->i * (addo * d->hamt - d->hamt * addo) + gamma * addo;
  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    if (key(i_nind)) {
      int mod = d->nmodmaxLabel(i_nind);
      for (int ii = 0; ii < d->nmodmax; ii++) {
        addo_fil += (parity_key(i_nind) ? 1.0 : -1.0) * d->i / d->coef_abs(i_nind) * ((double)parity_p_n * d->coef_lft[i_nind](mod, ii) * addo * d->qmdtc[i_nind] + d->coef_rht[i_nind](mod, ii) * d->qmdtc[i_nind] * addo);
      }
    } else {
      addo_fil += (parity_key(i_nind) ? -1.0 : 1.0) * d->i * d->coef_abs(i_nind) * ((double)parity_p_n * addo * d->qmdta[i_nind] - d->qmdta[i_nind] * addo);
    }
  }

  if (addo_fil.squaredNorm() > d->ferr || addo1.squaredNorm() > d->ferr) {
    int nddo_tmp;
#pragma omp critical
    {
      nddo_tmp = d->nddo;
      d->nddo++;
    }
    llint hash_val = generate_hash_value(key, d);
    tree->insert(hash_val, nddo_tmp);
    ddos1[nddo_tmp] = addo;
    ddos3[nddo_tmp] = addo1;
    keys1[nddo_tmp] = key;
  }
}

void filter_zero_p(DEOM_DATA *d, const nNNmat &ddos, const array_key_vec &keys, nNNmat &ddos1, array_key_vec &keys1, Trie *tree, int iado) {
  MatrixNcd addo;
  addo.resize(NSYS, NSYS);
  key_vec key(d->nind);
  key = keys[iado];
  addo = ddos[iado];

  if (is_valid(addo, 1e-15)) {
    int nddo_tmp;
#pragma omp critical
    {
      nddo_tmp = d->nddo;
      d->nddo++;
    }
    llint hash_val = generate_hash_value(key, d);
    tree->insert(hash_val, nddo_tmp);
    ddos1[nddo_tmp] = addo;
    keys1[nddo_tmp] = key;
  }
}