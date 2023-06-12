#pragma once
#include "deom.hpp"

int generate_key_plus(const key_vec &key0, key_vec &key1, const int pos,
                      const DEOM_DATA *d) {
  int sum = 0, hash_value = 1;
  key1 = key0;
  key1(pos) = true;
  for (int i = 0; i < d->nind; i++) {
    if (key1(i)) {
      sum += 1;
      hash_value += d->comb_list(i, sum) + d->comb_list(d->nind, sum - 1);
    }
  }
  return hash_value;
}

int generate_key_minus(const key_vec &key0, key_vec &key1, const int pos,
                       const DEOM_DATA *d) {
  int sum = 0, hash_value = 1;
  key1 = key0;
  key1(pos) = false;
  for (int i = 0; i < d->nind; i++) {
    if (key1(i)) {
      sum += 1;
      hash_value += d->comb_list(i, sum) + d->comb_list(d->nind, sum - 1);
    }
  }
  return hash_value;
}

int generate_hash_value_plus(key_vec &key0, const int pos,
                             const DEOM_DATA *d) {
  int sum = 0, hash_value = 1;
  key0(pos) = true;
  for (int i = 0; i < d->nind; i++) {
    if (key0(i)) {
      sum += 1;
      hash_value += d->comb_list(i, sum) + d->comb_list(d->nind, sum - 1);
    }
  }
  key0(pos) = false;
  return hash_value;
}

int generate_hash_value_minus(key_vec &key0, const int pos,
                              const DEOM_DATA *d) {
  int sum = 0, hash_value = 1;
  key0(pos) = false;
  for (int i = 0; i < d->nind; i++) {
    if (key0(i)) {
      sum += 1;
      hash_value += d->comb_list(i, sum) + d->comb_list(d->nind, sum - 1);
    }
  }
  key0(pos) = true;
  return hash_value;
}

int tier(const key_vec &key, const DEOM_DATA *d) {
  int tier = 0;
  for (int i = 0; i < d->nind; i++) tier += key(i);
  return tier;
}

llint gen_hash_value(const key_vec &key, const DEOM_DATA *d) {
  int sum = 0, hash_value = 1;
  for (int i = 0; i < d->nind; i++) {
    if (key(i)) {
      sum += 1;
      hash_value += d->comb_list(i, sum) + d->comb_list(d->nind, sum - 1);
    }
  }
  return (llint)hash_value;
}

void gen_parity(const DEOM_DATA *d, const key_vec &key,
                key_vec &parity_key, double *parity_p_n,
                complex<double> *gams) {
  bool parity = true;
  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    if (key(i_nind)) {
      *gams += d->expn(i_nind);
      parity ^= true;
    }
    parity_key(i_nind) = parity;
  }
  *parity_p_n = parity ? 1.0 : -1.0;
}

void gen_parity(const DEOM_DATA *d, const key_vec &key,
                key_vec &parity_key, double *parity_p_n) {
  bool parity = true;
  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    if (key(i_nind)) {
      parity ^= true;
    }
    parity_key(i_nind) = parity;
  }
  *parity_p_n = parity ? 1.0 : -1.0;
}

void gen_parity(const DEOM_DATA *d, const key_vec &key,
                double *parity_p_n) {
  bool parity = true;
  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    if (key(i_nind)) {
      parity ^= true;
    }
  }
  *parity_p_n = parity ? 1.0 : -1.0;
}