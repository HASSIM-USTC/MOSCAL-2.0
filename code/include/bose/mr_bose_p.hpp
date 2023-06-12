#pragma once
#include "algebra.hpp"
#include "deom.hpp"
#include "linearalgebra.hpp"

void construct_Mapping_p(DEOM_DATA *d, int iado) {
  key_vec key(d->keys[iado]), key_swell(d->nind);
  int hash_val = gen_hash_value(key, d);

  if (key.sum() < d->lmax) {
    for (int mp = 0; mp < d->nind; mp++) {
      hash_val = generate_key_plus(key, key_swell, mp, d) - 1;
      d->keys[hash_val] = key_swell;
    }
  }

  for (int mp = 0; mp < d->nind; mp++) {
    if (key(mp) > 0) {
      hash_val = generate_key_minus(key, key_swell, mp, d) - 1;
      d->keys[hash_val] = key_swell;
    }
  }
}

void rem_oprt_p(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
                const int iado, const char lcr) {}

void rem_cal_p(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
               const int iado) {
  int pos, tier = 0;
  complex<double> gams = 0.0;
  key_vec key(d->keys[iado]), key_swell(d->nind);

  for (int i = 0; i < d->nind; i++) {
    gams += (double)(key(i)) * d->expn(i);
    tier += key(i);
  }

  ddos1[iado] = -gams * total[iado];
  ddos1[iado] -= d->i * (d->hamt * total[iado] - total[iado] * d->hamt);

  for (int mp = 0; mp < d->nind; mp++) {
    double n = key(mp);
    if (key(mp) > 0) {
      pos = generate_key_minus(key, key_swell, mp, d) - 1;
      ddos1[iado] -=
          d->i * sqrt(n) *
          (d->qmdtc_l[mp] * total[pos] - total[pos] * d->qmdtc_r[mp]);
    }

    if (tier < d->lmax) {
      pos = generate_key_plus(key, key_swell, mp, d) - 1;
      ddos1[iado] -=
          d->i * sqrt(n + 1) *
          (d->qmdta_l[mp] * total[pos] - total[pos] * d->qmdta_r[mp]);
    }
  }
}

void rem_cal_1_p(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
                 const int iado) {}

void syl_cal_p(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
               const DEOM_SYL *syl, const int iado) {}

void syl_cal_1_p(nNNmat &ddos1, const nNNmat &total, const nNNmat &ddos2,
                 const DEOM_DATA *d, const DEOM_SYL *syl, const int iado) {}

void syl_cal_ddos_p(nNNmat &ddos1, const nNNmat &total, const nNNmat &ddos2,
                    const DEOM_DATA *d, const DEOM_SYL *syl, const int iado) {}

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