#pragma once

#include "algebra.hpp"
#include "deom.hpp"
#include "linearalgebra.hpp"

void construct_Mapping_p(DEOM_DATA *d, int iado) {
  key_vec key(d->keys[iado]), key_swell(d->nind);
  llint hash_val;
  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    if (!key(i_nind)) {
      if (tier(key, d) < d->lmax) {
        hash_val = generate_key_plus(key, key_swell, i_nind, d) - 1;
        if (hash_val != -1) d->keys[hash_val] = key_swell;
      }
    }
  }
}

void rem_oprt_p(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
                const int iado, const char lcr) {
  key_vec key(d->keys[iado]), key_swell(d->nind), parity_key(d->nind);
  double parity_p_n = 0;
  gen_parity(d, key, parity_key, &parity_p_n);
  ArrayXi g_index_list(d->nind);

  ddos1[iado].setZero();

  for (int mp = 0; mp < d->nmod; mp++) {
    if (lcr == 'l')
      ddos1[iado] += parity_p_n * (d->dipole.sdip[mp] * total[iado]);
    else if (lcr == 'r')
      ddos1[iado] += (total[iado] * d->dipole.sdip[mp]);
    else if (lcr == 'c')
      ddos1[iado] += (parity_p_n * d->dipole.sdip[mp] * total[iado] +
                      total[iado] * d->dipole.sdip[mp]);
  }

  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    if (key(i_nind)) {
      int pos = generate_key_minus(key, key_swell, i_nind, d) - 1;
      if (lcr == 'l') {
        ddos1[iado] -= (parity_key(i_nind) ? 1.0 : -1.0) * parity_p_n *
                       d->dipole.qmdtc_l[i_nind] * total[pos];
      } else if (lcr == 'r') {
        ddos1[iado] -= -(parity_key(i_nind) ? 1.0 : -1.0) * total[pos] *
                       d->dipole.qmdtc_r[i_nind];
      } else if (lcr == 'c') {
        ddos1[iado] -=
            (parity_key(i_nind) ? 1.0 : -1.0) *
            (parity_p_n * d->dipole.qmdtc_l[i_nind] * total[pos] -
             total[pos] * d->dipole.qmdtc_r[i_nind]);
      }
    } else if (tier(key, d) < d->lmax) {
      int pos = generate_key_plus(key, key_swell, i_nind, d) - 1;
      if (lcr == 'l') {
        ddos1[iado] -=
            (parity_key(i_nind) ? 1.0 : -1.0) *
            (parity_p_n * d->dipole.qmdta_l[i_nind] * total[pos]);
      } else if (lcr == 'r') {
        ddos1[iado] -= (parity_key(i_nind) ? 1.0 : -1.0) * total[pos] *
                       d->dipole.qmdta_r[i_nind];
      } else if (lcr == 'c') {
        ddos1[iado] -=
            (parity_key(i_nind) ? 1.0 : -1.0) *
            (parity_p_n * d->dipole.qmdta_l[i_nind] * total[pos] +
             total[pos] * d->dipole.qmdta_r[i_nind]);
      }
    }
  }

  init_nNNmat(ddos1[iado]);
}

void rem_cal_p(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
               const int iado) {
  key_vec key(d->keys[iado]), key_swell(d->nind), parity_key(d->nind);
  int pos, pos1;
  complex<double> gams = 0;
  double parity_p_n = 0;
  gen_parity(d, key, parity_key, &parity_p_n, &gams);

  ddos1[iado] = -(gams * total[iado]);
  ddos1[iado] -= d->i * (d->hamt * total[iado] - total[iado] * d->hamt);

  for (int mp = 0; mp < d->nind; mp++) {
    if (key(mp)) {
      pos = generate_key_minus(key, key_swell, mp, d) - 1;
      ddos1[iado] -= (double)(parity_key(mp) ? 1 : -1) * d->i * d->alp1 *
                     ((double)parity_p_n * d->qmdtc_l[mp] * total[pos] +
                      total[pos] * d->qmdtc_r[mp]);

      for (int np = 0; np < d->nind; np++) {
        if (key_swell(np)) {
          if (d->qmdt2c_equal_0(mp, np)) {
            pos1 = generate_hash_value_minus(key_swell, np, d) - 1;
            if (np < mp) {
              ddos1[iado] -=
                  (double)(parity_key(mp) ? 1 : -1) *
                  (double)(parity_key(np) ? 1 : -1) * d->i *
                  (d->qmdt2c_l[INDEX2(mp, np, d->nind)] * total[pos1] -
                   total[pos1] * d->qmdt2c_r[INDEX2(mp, np, d->nind)]);
            }
          }
        } else {
          if (d->qmdt2b_equal_0(mp, np)) {
            pos1 = generate_hash_value_plus(key_swell, np, d) - 1;
            if (np < mp) {
              ddos1[iado] -=
                  (double)(parity_key(mp) ? 1 : -1) *
                  (double)(parity_key(np) ? 1 : -1) * d->i *
                  (d->qmdt2b_l[INDEX2(mp, np, d->nind)] * total[pos1] +
                   total[pos1] * d->qmdt2b_r[INDEX2(mp, np, d->nind)]);
            } else {
              ddos1[iado] +=
                  (double)(parity_key(mp) ? 1 : -1) *
                  (double)(parity_key(np) ? 1 : -1) * d->i *
                  (d->qmdt2b_l[INDEX2(mp, np, d->nind)] * total[pos1] +
                   total[pos1] * d->qmdt2b_r[INDEX2(mp, np, d->nind)]);
            }
          }
        }
      }
    } else if (tier(key, d) < d->lmax) {
      pos = generate_key_plus(key, key_swell, mp, d) - 1;
      ddos1[iado] -= (double)(parity_key(mp) ? 1 : -1) * d->i * d->alp1 *
                     ((double)parity_p_n * d->qmdta_l[mp] * total[pos] -
                      total[pos] * d->qmdta_r[mp]);
      if (tier(key, d) + 1 < d->lmax) {
        for (int np = 0; np < d->nind; np++) {
          if (d->qmdt2a_equal_0(mp, np)) {
            if (!key_swell(np)) {
              pos1 = generate_hash_value_plus(key_swell, np, d) - 1;
              if (np < mp) {
                ddos1[iado] -=
                    0.5 * (double)(parity_key(mp) ? 1 : -1) *
                    (double)(parity_key(np) ? 1 : -1) * d->i *
                    (d->qmdt2a_l[INDEX2(mp, np, d->nind)] * total[pos1] -
                     total[pos1] * d->qmdt2a_r[INDEX2(mp, np, d->nind)]);
              } else if (np > mp) {
                ddos1[iado] +=
                    0.5 * (double)(parity_key(mp) ? 1 : -1) *
                    (double)(parity_key(np) ? 1 : -1) * d->i *
                    (d->qmdt2a_l[INDEX2(mp, np, d->nind)] * total[pos1] -
                     total[pos1] * d->qmdt2a_r[INDEX2(mp, np, d->nind)]);
              }
            }
          }
        }
      }
    }
  }

  init_nNNmat(ddos1[iado]);
}

void rem_cal_1_p(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
                 const int iado) {}

void syl_cal_p(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
               const DEOM_SYL *syl, const int iado) {
  key_vec key(d->keys[iado]), key_swell(d->nind), parity_key(d->nind);
  int pos, pos1;
  ArrayXi key_syl(syl->nind);
  double parity_p_n = 0;
  gen_syl_key(key_syl, key, syl);
  gen_parity(d, key, parity_key, &parity_p_n);

  ddos1[iado] = syl->OMG * total[iado];

  for (int mp = 0; mp < d->nind; mp++) {
    if (key(mp)) {
      pos = generate_key_minus(key, key_swell, mp, d) - 1;
      ddos1[iado] -= (double)(parity_key(mp) ? 1 : -1) * d->i * d->alp1 *
                     ((double)parity_p_n * d->qmdtc_l[mp] * total[pos] +
                      total[pos] * d->qmdtc_r[mp]);

      for (int np = 0; np < d->nind; np++) {
        if (key_swell(np)) {
          if (d->qmdt2c_equal_0(mp, np)) {
            pos1 = generate_hash_value_minus(key_swell, np, d) - 1;
            if (np < mp)
              ddos1[iado] -=
                  (double)(parity_key(mp) ? 1 : -1) *
                  (double)(parity_key(np) ? 1 : -1) * d->i *
                  (d->qmdt2c_l[INDEX2(mp, np, d->nind)] * total[pos1] -
                   total[pos1] * d->qmdt2c_r[INDEX2(mp, np, d->nind)]);
          }
        } else {
          if (d->qmdt2b_equal_0(mp, np)) {
            pos1 = generate_hash_value_plus(key_swell, np, d) - 1;
            if (np < mp)
              ddos1[iado] -=
                  (double)(parity_key(mp) ? 1 : -1) *
                  (double)(parity_key(np) ? 1 : -1) * d->i *
                  (d->qmdt2b_l[INDEX2(mp, np, d->nind)] * total[pos1] +
                   total[pos1] * d->qmdt2b_r[INDEX2(mp, np, d->nind)]);
            else
              ddos1[iado] +=
                  (double)(parity_key(mp) ? 1 : -1) *
                  (double)(parity_key(np) ? 1 : -1) * d->i *
                  (d->qmdt2b_l[INDEX2(mp, np, d->nind)] * total[pos1] +
                   total[pos1] * d->qmdt2b_r[INDEX2(mp, np, d->nind)]);
          }
        }
      }
    } else if (tier(key, d) < d->lmax) {
      pos = generate_key_plus(key, key_swell, mp, d) - 1;
      ddos1[iado] -= (double)(parity_key(mp) ? 1 : -1) * d->i * d->alp1 *
                     ((double)parity_p_n * d->qmdta_l[mp] * total[pos] -
                      total[pos] * d->qmdta_r[mp]);
      if (tier(key, d) + 1 < d->lmax) {
        for (int np = 0; np < d->nind; np++) {
          if (d->qmdt2a_equal_0(mp, np)) {
            if (!key_swell(np)) {
              pos1 = generate_hash_value_plus(key_swell, np, d) - 1;
              if (np < mp)
                ddos1[iado] -=
                    0.5 * (double)(parity_key(mp) ? 1 : -1) *
                    (double)(parity_key(np) ? 1 : -1) * d->i *
                    (d->qmdt2a_l[INDEX2(mp, np, d->nind)] * total[pos1] -
                     total[pos1] * d->qmdt2a_r[INDEX2(mp, np, d->nind)]);
              else if (np > mp)
                ddos1[iado] +=
                    0.5 * (double)(parity_key(mp) ? 1 : -1) *
                    (double)(parity_key(np) ? 1 : -1) * d->i *
                    (d->qmdt2a_l[INDEX2(mp, np, d->nind)] * total[pos1] -
                     total[pos1] * d->qmdt2a_r[INDEX2(mp, np, d->nind)]);
            }
          }
        }
      }
    }
  }

  init_nNNmat(ddos1[iado]);

  pos = generate_hash_value_syl(key_syl, d, syl) - 1;
  syl_complex_storge(ddos1[iado], syl, pos, d->ferr, ddos1[iado]);
}

void syl_cal_1_p(nNNmat &ddos1, const nNNmat &total, const nNNmat &ddos2,
                 const DEOM_DATA *d, const DEOM_SYL *syl, const int iado) {
  key_vec key(d->keys[iado]), key_swell(d->nind), parity_key(d->nind);
  int pos, pos1;
  double parity_p_n = 0;
  ArrayXi key_syl(syl->nind);
  gen_syl_key(key_syl, key, syl);
  gen_parity(d, key, parity_key, &parity_p_n);

  ddos1[iado] = syl->OMG * total[iado] + ddos2[iado];

  for (int mp = 0; mp < d->nind; mp++) {
    if (key(mp)) {
      pos = generate_key_minus(key, key_swell, mp, d) - 1;
      ddos1[iado] -= (double)(parity_key(mp) ? 1 : -1) * d->i * d->alp1 *
                     ((double)parity_p_n * d->qmdtc_l[mp] * total[pos] +
                      total[pos] * d->qmdtc_r[mp]);

      for (int np = 0; np < d->nind; np++) {
        if (key_swell(np)) {
          if (d->qmdt2c_equal_0(mp, np)) {
            pos1 = generate_hash_value_minus(key_swell, np, d) - 1;
            if (np < mp)
              ddos1[iado] -=
                  (double)(parity_key(mp) ? 1 : -1) *
                  (double)(parity_key(np) ? 1 : -1) * d->i *
                  (d->qmdt2c_l[INDEX2(mp, np, d->nind)] * total[pos1] -
                   total[pos1] * d->qmdt2c_r[INDEX2(mp, np, d->nind)]);
          }
        } else {
          if (d->qmdt2b_equal_0(mp, np)) {
            pos1 = generate_hash_value_plus(key_swell, np, d) - 1;
            if (np < mp)
              ddos1[iado] -=
                  (double)(parity_key(mp) ? 1 : -1) *
                  (double)(parity_key(np) ? 1 : -1) * d->i *
                  (d->qmdt2b_l[INDEX2(mp, np, d->nind)] * total[pos1] +
                   total[pos1] * d->qmdt2b_r[INDEX2(mp, np, d->nind)]);
            else
              ddos1[iado] +=
                  (double)(parity_key(mp) ? 1 : -1) *
                  (double)(parity_key(np) ? 1 : -1) * d->i *
                  (d->qmdt2b_l[INDEX2(mp, np, d->nind)] * total[pos1] +
                   total[pos1] * d->qmdt2b_r[INDEX2(mp, np, d->nind)]);
          }
        }
      }
    } else if (tier(key, d) < d->lmax) {
      pos = generate_key_plus(key, key_swell, mp, d) - 1;
      ddos1[iado] -= (double)(parity_key(mp) ? 1 : -1) * d->i * d->alp1 *
                     ((double)parity_p_n * d->qmdta_l[mp] * total[pos] -
                      total[pos] * d->qmdta_r[mp]);
      if (tier(key, d) + 1 < d->lmax) {
        for (int np = 0; np < d->nind; np++) {
          if (d->qmdt2a_equal_0(mp, np)) {
            if (!key_swell(np)) {
              pos1 = generate_hash_value_plus(key_swell, np, d) - 1;
              if (np < mp)
                ddos1[iado] -=
                    0.5 * (double)(parity_key(mp) ? 1 : -1) *
                    (double)(parity_key(np) ? 1 : -1) * d->i *
                    (d->qmdt2a_l[INDEX2(mp, np, d->nind)] * total[pos1] -
                     total[pos1] * d->qmdt2a_r[INDEX2(mp, np, d->nind)]);
              else if (np > mp)
                ddos1[iado] +=
                    0.5 * (double)(parity_key(mp) ? 1 : -1) *
                    (double)(parity_key(np) ? 1 : -1) * d->i *
                    (d->qmdt2a_l[INDEX2(mp, np, d->nind)] * total[pos1] -
                     total[pos1] * d->qmdt2a_r[INDEX2(mp, np, d->nind)]);
            }
          }
        }
      }
    }
  }

  init_nNNmat(ddos1[iado]);

  pos = generate_hash_value_syl(key_syl, d, syl) - 1;
  syl_complex_storge(ddos1[iado], syl, pos, d->ferr, ddos1[iado]);
}

void syl_cal_ddos_p(nNNmat &ddos1, const nNNmat &total,
                    const nNNmat &ddos2, const DEOM_DATA *d,
                    const DEOM_SYL *syl, const int iado) {}

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
