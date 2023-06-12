#pragma once

#include "algebra.hpp"
#include "deom.hpp"

void construct_Mapping_filter_p_no_add(DEOM_DATA *d, Trie *tree, int iado) {
  key_vec key(d->keys[iado]), key_swell(d->nind);
  llint hash_val = gen_hash_value(key, d);
  int pos_get = tree->find(hash_val);
  d->g_index_list_pos[pos_get](0) = pos_get + 1;

  for (int mp = 0; mp < d->nind; mp++) {
    if (key.sum() < d->lmax) {
      hash_val = generate_key_plus(key, key_swell, mp, d);
      int pos = tree->find(hash_val);
      if (pos != -1)
        d->g_index_list_pos[pos](mp + 1) = pos_get + 1;
    }

    if (key(mp) > 0) {
      hash_val = generate_key_minus(key, key_swell, mp, d);
      int pos = tree->find(hash_val);
      if (pos != -1)
        d->g_index_list_pos[pos](mp + d->nind + 1) = pos_get + 1;
    }
  }
}

void construct_Mapping_filter_p(nNNmat &ddos, DEOM_DATA *d, Trie *tree, int iado) {
}

void filter_p(DEOM_DATA *d, const nNNmat &ddos, const nXvecb &keys, nNNmat &ddos1, nXvecb &keys1, Trie *tree, int iado) {
}

void filter_2_p(DEOM_DATA *d, const nNNmat &ddos, const nNNmat &ddos2, const nXvecb &keys, nNNmat &ddos1, nNNmat &ddos3, nXvecb &keys1, Trie *tree, int iado) {
}

void filter_hei_p(DEOM_DATA *d, const nNNmat &ddos, const nNNmat &ddos2, const nXvecb &keys, nNNmat &ddos1, nNNmat &ddos3, nXvecb &keys1, Trie *tree, int iado) {
}

void filter_zero_p(DEOM_DATA *d, const nNNmat &ddos, const nXvecb &keys, nNNmat &ddos1, nXvecb &keys1, Trie *tree, int iado) {
}