#pragma once
#include "deom.hpp"

int generate_key_plus(const key_vec &key0, key_vec &key1, const int pos, const DEOM_DATA *d) {
  int sum = 0, hash_value = 1;
  key1 = key0;
  key1(pos) += 1;
  for (int i = 0; i < d->nind; i++) {
    sum += key1(i);
    hash_value += d->comb_list(sum + i, i + 1);
  }
  return hash_value;
}

int generate_key_minus(const key_vec &key0, key_vec &key1, const int pos, const DEOM_DATA *d) {
  int sum = 0, hash_value = 1;
  key1 = key0;
  key1(pos) += -1;
  for (int i = 0; i < d->nind; i++) {
    sum += key1(i);
    hash_value += d->comb_list(sum + i, i + 1);
  }
  return hash_value;
}

int generate_hash_value_plus(key_vec &key0, const int pos, const DEOM_DATA *d) {
  int sum = 0, hash_value = 1;
  key0(pos) += 1;
  for (int i = 0; i < d->nind; i++) {
    sum += key0(i);
    hash_value += d->comb_list(sum + i, i + 1);
  }
  key0(pos) -= 1;
  return hash_value;
}

int generate_hash_value_minus(key_vec &key0, const int pos, const DEOM_DATA *d) {
  int sum = 0, hash_value = 1;
  key0(pos) -= 1;
  for (int i = 0; i < d->nind; i++) {
    sum += key0(i);
    hash_value += d->comb_list(sum + i, i + 1);
  }
  key0(pos) += 1;
  return hash_value;
}

int tier(const key_vec &key, const DEOM_DATA *d) {
  int tier = 0;
  for (int i = 0; i < d->nind; i++)
    tier += key(i);
  return tier;
}

llint gen_hash_value(const key_vec &key, const DEOM_DATA *d) {
  int sum = 0, hash_value = 1;
  for (int i = 0; i < d->nind; i++) {
    sum += key(i);
    hash_value += d->comb_list(sum + i, i + 1);
  }
  return (llint)hash_value;
}
