#pragma once

#include "algebra.hpp"
#include "deom.hpp"

void construct_Mapping_filter_p_no_add(DEOM_DATA *d, Trie *tree, int iado) {
  key_vec key(d->keys[iado]), key_swell(d->nind), key_swell1(d->nind);
  llint hash_val = gen_hash_value(key, d);
  int pos, pos_get = tree->find(hash_val);
  d->g_index_list_pos[pos_get](0) = pos_get + 1;

  for (int mp1 = 0; mp1 < d->nind; mp1++) {
    hash_val = generate_key_plus(key, key_swell, mp1, d);
    if (key.sum() < d->lmax) {
      pos = tree->find(hash_val);
      if (pos != -1)
        d->g_index_list_pos[pos](mp1 + 1) = pos_get + 1;
    }

    if (key.sum() < d->lmax - 1) {
      for (int mp2 = 0; mp2 < d->nind; mp2++) {
        hash_val = generate_key_plus(key_swell, key_swell1, mp2, d);
        int pos = tree->find(hash_val);
        if (pos != -1)
          d->g_index_list_pos[pos](2 * d->nind + 1 + mp1 * d->nind + mp2) = pos_get + 1;
      }
    }

    for (int mp2 = 0; mp2 < d->nind; mp2++) {
      if (key_swell(mp2) > 0) {
        hash_val = generate_key_minus(key_swell, key_swell1, mp2, d);
        int pos = tree->find(hash_val);
        if (pos != -1)
          d->g_index_list_pos[pos](2 * d->nind + 1 + d->nind * d->nind + mp1 * d->nind + mp2) = pos_get + 1;
      }
    }
  }

  for (int mp1 = 0; mp1 < d->nind; mp1++) {
    if (key(mp1) > 0) {
      hash_val = generate_key_minus(key, key_swell, mp1, d);
      int pos = tree->find(hash_val);
      if (pos != -1)
        d->g_index_list_pos[pos](mp1 + d->nind + 1) = pos_get + 1;

      for (int mp2 = 0; mp2 < d->nind; mp2++) {
        if (key_swell(mp2) > 0) {
          hash_val = generate_key_minus(key_swell, key_swell1, mp2, d);
          int pos = tree->find(hash_val);
          if (pos != -1)
            d->g_index_list_pos[pos](2 * d->nind + 1 + 2 * d->nind * d->nind + mp1 * d->nind + mp2) = pos_get + 1;
        }
      }
    }
  }
}

void construct_Mapping_filter_p(nNNmat &ddos, DEOM_DATA *d, Trie *tree, int iado) {
  key_vec key(d->nind), key_swell(d->nind), key_swell1(d->nind);
  llint hash_val;
  int pos, pos_get;

  key.noalias() = d->keys[iado];
  hash_val = generate_hash_value(key, d);
  if (is_valid(ddos, iado, d->ferr) || hash_val == 1) {
    pos_get = tree->find(hash_val);
    d->g_index_list_pos[pos_get](0) = pos_get + 1;

    for (int mp1 = 0; mp1 < d->nind; mp1++) {
      if (key.sum() < d->lmax) {
        hash_val = generate_key_plus(key, key_swell, mp1, d);
        pos = tree->find(hash_val);
        if (pos == -1) {
          int nddo_tmp;
#pragma omp critical
          {
            nddo_tmp = d->nddo;
            d->nddo++;
          }
          const auto [rank_find, success] = tree->try_insert(hash_val, nddo_tmp);
          if (success) {
            d->g_index_list_pos[rank_find](mp1 + 1) = pos_get + 1;
            d->keys[rank_find] = key_swell;
          } else {
            pos = tree->find(hash_val);
            d->g_index_list_pos[pos](mp1 + 1) = pos_get + 1;
            d->keys[nddo_tmp].setZero();
          }
        } else
          d->g_index_list_pos[pos](mp1 + 1) = pos_get + 1;

        for (int mp2 = 0; mp2 < d->nind; mp2++) {
          if (key.sum() < d->lmax - 1) {
            hash_val = generate_key_plus(key_swell, key_swell1, mp2, d);
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
                d->g_index_list_pos[rank_find](2 * d->nind + 1 + mp1 * d->nind + mp2) = pos_get + 1;
                d->keys[rank_find] = key_swell1;
              } else {
                pos = tree->find(hash_val);
                d->g_index_list_pos[pos](2 * d->nind + 1 + mp1 * d->nind + mp2) = pos_get + 1;
                d->keys[nddo_tmp].setZero();
              }
            } else
              d->g_index_list_pos[pos](2 * d->nind + 1 + mp1 * d->nind + mp2) = pos_get + 1;
          }

          if (key_swell(mp2) > 0) {
            hash_val = generate_key_minus(key_swell, key_swell1, mp2, d);
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
                d->g_index_list_pos[rank_find](2 * d->nind + 1 + d->nind * d->nind + mp1 * d->nind + mp2) = pos_get + 1;
                d->keys[rank_find] = key_swell1;
              } else {
                pos = tree->find(hash_val);
                d->g_index_list_pos[pos](2 * d->nind + 1 + d->nind * d->nind + mp1 * d->nind + mp2) = pos_get + 1;
                d->keys[nddo_tmp].setZero();
              }
            } else
              d->g_index_list_pos[pos](2 * d->nind + 1 + d->nind * d->nind + mp1 * d->nind + mp2) = pos_get + 1;
          }
        }
      }

      if (key(mp1) > 0) {
        hash_val = generate_key_minus(key, key_swell, mp1, d);
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
            d->g_index_list_pos[rank_find](mp1 + d->nind + 1) = pos_get + 1;
            d->keys[rank_find] = key_swell;
          } else {
            pos = tree->find(hash_val);
            d->g_index_list_pos[pos](mp1 + d->nind + 1) = pos_get + 1;
            d->keys[nddo_tmp].setZero();
          }
        } else
          d->g_index_list_pos[pos](mp1 + d->nind + 1) = pos_get + 1;

        for (int mp2 = 0; mp2 < d->nind; mp2++) {
          if (key_swell(mp2) > 0) {
            hash_val = generate_key_minus(key_swell, key_swell1, mp2, d);
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
                d->g_index_list_pos[rank_find](2 * d->nind + 1 + 2 * d->nind * d->nind + mp1 * d->nind + mp2) = pos_get + 1;
                d->keys[rank_find] = key_swell1;
              } else {
                pos = tree->find(hash_val);
                d->g_index_list_pos[pos](2 * d->nind + 1 + 2 * d->nind * d->nind + mp1 * d->nind + mp2) = pos_get + 1;
                d->keys[nddo_tmp].setZero();
              }
            } else
              d->g_index_list_pos[pos](2 * d->nind + 1 + 2 * d->nind * d->nind + mp1 * d->nind + mp2) = pos_get + 1;
          }
        }
      }
    }
  }
}

void filter_p(DEOM_DATA *d, const nNNmat &ddos, const array_key_vec &keys, nNNmat &ddos1, array_key_vec &keys1, Trie *tree, int iado) {
  if (is_valid(ddos[iado], d->ferr)) {
    int nddo_tmp;
#pragma omp critical
    {
      nddo_tmp = d->nddo;
      d->nddo++;
    }
    key_vec key(d->nind);
    key = keys[iado];
    llint hash_val = generate_hash_value(key, d);
    tree->insert(hash_val, nddo_tmp);
    ddos1[nddo_tmp] = ddos[iado];
    keys1[nddo_tmp] = key;
  }
}

void filter_2_p(DEOM_DATA *d, const nNNmat &ddos, const nNNmat &ddos2, const array_key_vec &keys, nNNmat &ddos1, nNNmat &ddos3, array_key_vec &keys1, Trie *tree, int iado) {
  key_vec key(d->nind);
  key = keys[iado];

  complex<double> absl = 2 * d->ham1_abs;
  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    if (key(i_nind)) {
      absl += abs(d->expn(i_nind));
      absl += 2 * d->qmd1_abs(i_nind) / d->coef_abs(i_nind);
    } else {
      absl += 2 * d->qmd1_abs(i_nind) * d->coef_abs(i_nind);
    }
  }

  if (is_valid(ddos[iado], absl, d->ferr)) {
    int nddo_tmp;
#pragma omp critical
    {
      nddo_tmp = d->nddo;
      d->nddo++;
    }
    llint hash_val = generate_hash_value(key, d);
    tree->insert(hash_val, nddo_tmp);
    ddos1[nddo_tmp] = ddos[iado];
    ddos3[nddo_tmp] = ddos2[iado];
    keys1[nddo_tmp] = key;
  }
}

void filter_hei_p(DEOM_DATA *d, const nNNmat &ddos, const nNNmat &ddos2, const array_key_vec &keys, nNNmat &ddos1, nNNmat &ddos3, array_key_vec &keys1, Trie *tree, int iado) {
  key_vec key(d->nind);
  key = keys[iado];

  complex<double> absl = 2 * d->ham1_abs;
  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    if (key(i_nind)) {
      absl += abs(d->expn(i_nind));
      absl += 2 * d->qmd1_abs(i_nind) / d->coef_abs(i_nind);
    } else {
      absl += 2 * d->qmd1_abs(i_nind) * d->coef_abs(i_nind);
    }
  }
  if (is_valid(ddos[iado], absl, d->ferr) || ddos2[iado].squaredNorm() > d->ferr) {
    int nddo_tmp;
#pragma omp critical
    {
      nddo_tmp = d->nddo;
      d->nddo++;
    }
    llint hash_val = generate_hash_value(key, d);
    tree->insert(hash_val, nddo_tmp);
    ddos1[nddo_tmp] = ddos[iado];
    ddos3[nddo_tmp] = ddos2[iado];
    keys1[nddo_tmp] = key;
  }
}

void filter_zero_p(DEOM_DATA *d, const nNNmat &ddos, const array_key_vec &keys, nNNmat &ddos1, array_key_vec &keys1, Trie *tree, int iado) {
  key_vec key(d->nind);
  key = keys[iado];

  if (is_valid(ddos[iado], 1e-15)) {
    int nddo_tmp;
#pragma omp critical
    {
      nddo_tmp = d->nddo;
      d->nddo++;
    }
    llint hash_val = generate_hash_value(key, d);
    tree->insert(hash_val, nddo_tmp);
    ddos1[nddo_tmp] = ddos[iado];
    keys1[nddo_tmp] = key;
  }
}