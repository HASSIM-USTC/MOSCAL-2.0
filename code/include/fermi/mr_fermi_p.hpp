#pragma once
#include "algebra.hpp"
#include "deom.hpp"
#include "linearalgebra.hpp"

void construct_Mapping_p(DEOM_DATA *d, int iado) {
  key_vec key(d->keys[iado]), key_swell(d->nind);
  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    if (!key(i_nind)) {
      if (tier(key, d) < d->lmax) {
        llint hash_val = generate_key_plus(key, key_swell, i_nind, d) - 1;
        if (hash_val != -1) {
          d->keys[hash_val] = key_swell;
        }
      }
    }
  }
}

complex<double> caculate_trace(DEOM_DATA *d) {
  MatrixXcd sum_Matrix = MatrixXcd::Zero(NSYS, NSYS);
  for (int nmod = 0; nmod < d->nmod; nmod++)
    sum_Matrix += d->dipole1.sdip[nmod] * d->ddos[0];

  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    int pos = generate_key_plus(d->zerokey, d->emptykey, i_nind, d) - 1;
    sum_Matrix += d->dipole1.qmdta_l[i_nind] * d->ddos[pos];
  }

  return sum_Matrix.trace();
}

void construct_Mapping_filter_p(DEOM_DATA *d, Trie *tree, int iado) {
  key_vec key(d->keys[iado]), key_swell(d->nind);
  llint hash_val = gen_hash_value(key, d);
  int pos_get = tree->find(hash_val);
  d->g_index_list_pos[pos_get](0) = pos_get + 1;

  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    if (key(i_nind)) {
      hash_val = generate_key_minus(key, key_swell, i_nind, d);
      set_ddos(key_swell, i_nind + 1, pos_get, hash_val, tree, d);
    } else {
      if (tier(key, d) < d->lmax) {
        hash_val = generate_key_plus(key, key_swell, i_nind, d);
        set_ddos(key_swell, i_nind + 1, pos_get, hash_val, tree, d);
      }
    }
  }
}

void construct_Mapping_filter_p(nNNmat &ddos, DEOM_DATA *d, Trie *tree,
                                int iado) {
  key_vec key(d->keys[iado]), key_swell(d->nind);
  llint hash_val = gen_hash_value(key, d);
  bool flag = is_valid(key, ddos[iado], d);

  if (flag || hash_val == 1) {
    int pos_get = tree->find(hash_val);
    d->g_index_list_pos[pos_get](0) = pos_get + 1;
    for (int i_nind = 0; i_nind < d->nind; i_nind++) {
      if (key(i_nind)) {
        hash_val = generate_key_minus(key, key_swell, i_nind, d);
        set_ddos(key_swell, i_nind + 1, pos_get, hash_val, tree, d);
      } else {
        if (tier(key, d) < d->lmax) {
          hash_val = generate_key_plus(key, key_swell, i_nind, d);
          set_ddos(key_swell, i_nind + 1, pos_get, hash_val, tree, d);
        }
      }
    }
  }
}

complex<double> caculate_trace(DEOM_DATA *d, Trie *tree) {
  MatrixXcd sum_Matrix = MatrixXcd::Zero(NSYS, NSYS);

  int pos = tree->find(1);
  if (pos != -1)
    for (int nmod = 0; nmod < d->nmod; nmod++)
      sum_Matrix += d->dipole1.sdip[nmod] * d->ddos[pos];

  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    llint hash_val = generate_key_plus(d->zerokey, d->emptykey, i_nind, d);
    pos = tree->find(hash_val);
    if (pos != -1) {
      sum_Matrix += d->dipole1.qmdta_l[i_nind] * d->ddos[pos];
    }
  }

  return sum_Matrix.trace();
}