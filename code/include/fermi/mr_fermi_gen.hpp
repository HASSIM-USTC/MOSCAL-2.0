#pragma once
#include "algebra.hpp"
#include "cfw.hpp"
#include "deom.hpp"
#include "linearalgebra.hpp"

// note in the noise spectrum calculation. Only one case, lcr = 'l', is
// used.
template <typename FUN>
void rem_oprt_gen(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
                  const int iado, const char lcr, FUN gen_g_index) {
  key_vec key(d->keys[iado]), parity_key(d->nind);
  double parity_p_n = 0;
  gen_parity(d, key, parity_key, &parity_p_n);
  ArrayXi g_index_list(d->nind + 1);
  gen_g_index(d, key, g_index_list, iado);

  ddos1[iado].setZero();

  for (int i_nind = 0; i_nind < d->nmod; i_nind++) {
    if (lcr == 'l')
      ddos1[iado] += parity_p_n * d->dipole.sdip[i_nind] * total[iado];
    else if (lcr == 'r')
      ddos1[iado] += total[iado] * d->dipole.sdip[i_nind];
    else if (lcr == 'c')
      ddos1[iado] += parity_p_n * d->dipole.sdip[i_nind] * total[iado] +
                     total[iado] * d->dipole.sdip[i_nind];
  }

  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    int pos = g_index_list(i_nind + 1) - 1;
    if (pos != -1) {
      if (key(i_nind)) {
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
      } else {
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
  }

  init_nNNmat(ddos1[iado]);
}

template <typename FUN>
void rem_cal_gen(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
                 const int iado, FUN gen_g_index) {
  key_vec key(d->keys[iado]), parity_key(d->nind);
  complex<double> gams = 0;
  double parity_p_n = 0;
  gen_parity(d, key, parity_key, &parity_p_n, &gams);
  ArrayXi g_index_list(d->nind + 1);
  gen_g_index(d, key, g_index_list, iado);

  ddos1[iado].setZero();

  int pos = g_index_list(0) - 1;
  if (pos != -1) {
    ddos1[iado] = -(gams * total[pos]);
    ddos1[iado] -= d->i * (d->hamt * total[pos] - total[pos] * d->hamt);
  }

  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    pos = g_index_list(i_nind + 1) - 1;
    if (pos != -1) {
      if (key(i_nind)) {
        ddos1[iado] -= (parity_key(i_nind) ? 1.0 : -1.0) * d->i *
                       (parity_p_n * d->qmdtc_l[i_nind] * total[pos] +
                        total[pos] * d->qmdtc_r[i_nind]);
      } else {
        ddos1[iado] -= (parity_key(i_nind) ? 1.0 : -1.0) * d->i *
                       (parity_p_n * d->qmdta_l[i_nind] * total[pos] -
                        total[pos] * d->qmdta_r[i_nind]);
      }
    }
  }

  init_nNNmat(ddos1[iado]);
}

template <typename FUN>
void rem_cal_hei_gen(nNNmat &ddos1, const nNNmat &total,
                     const DEOM_DATA *d, const int iado, FUN gen_g_index) {
  key_vec key(d->keys[iado]), parity_key(d->nind);
  complex<double> gams = 0;
  double parity_p_n = 0;
  gen_parity(d, key, parity_key, &parity_p_n, &gams);
  ArrayXi g_index_list(d->nind + 1);
  gen_g_index(d, key, g_index_list, iado);

  ddos1[iado].setZero();

  int pos = g_index_list(0) - 1;
  if (pos != -1) {
    ddos1[iado] = -gams * total[pos];
    ddos1[iado] -= d->i * (total[pos] * d->hamt - d->hamt * total[pos]);
  }

  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    pos = g_index_list(i_nind + 1) - 1;
    if (pos != -1) {
      if (key(i_nind)) {
        ddos1[iado] -= (parity_key(i_nind) ? 1.0 : -1.0) * d->i *
                       (parity_p_n * total[pos] * d->qmdta_l[i_nind] +
                        d->qmdta_r[i_nind] * total[pos]);
      } else {
        ddos1[iado] -= (parity_key(i_nind) ? 1.0 : -1.0) * d->i *
                       (parity_p_n * total[pos] * d->qmdtc_l[i_nind] -
                        d->qmdtc_r[i_nind] * total[pos]);
      }
    }
  }

  init_nNNmat(ddos1[iado]);
}

template <typename FUN>
void rem_cal_1_gen(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
                   const int iado, FUN gen_g_index) {
  key_vec key(d->keys[iado]), parity_key(d->nind);
  complex<double> gams = 0;
  double parity_p_n = 0;
  gen_parity(d, key, parity_key, &parity_p_n, &gams);
  ArrayXi g_index_list(d->nind + 1);
  gen_g_index(d, key, g_index_list, iado);

  ddos1[iado].setZero();

  int pos = g_index_list(0) - 1;
  if (pos != -1) {
    ddos1[iado] = -gams * total[pos];
    ddos1[iado] -= d->i * (d->hamt * total[pos] - total[pos] * d->hamt);
  }

  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    pos = g_index_list(i_nind + 1) - 1;
    if (pos != -1) {
      if (key(i_nind)) {
        ddos1[iado] -= (parity_key(i_nind) ? 1.0 : -1.0) * d->i *
                       (parity_p_n * d->qmdtc_l[i_nind] * total[pos] -
                        total[pos] * d->qmdtc_r[i_nind]);
      } else {
        ddos1[iado] -= (parity_key(i_nind) ? 1.0 : -1.0) * d->i *
                       (parity_p_n * d->qmdta_l[i_nind] * total[pos] +
                        total[pos] * d->qmdta_r[i_nind]);
      }
    }
  }

  init_nNNmat(ddos1[iado]);
}

template <typename FUN>
void rem_cal_1_hei_gen(nNNmat &ddos1, const nNNmat &total,
                       const DEOM_DATA *d, const int iado,
                       FUN gen_g_index) {
  key_vec key(d->keys[iado]), parity_key(d->nind);
  complex<double> gams = 0;
  double parity_p_n = 0;
  gen_parity(d, key, parity_key, &parity_p_n, &gams);
  ArrayXi g_index_list(d->nind + 1);
  gen_g_index(d, key, g_index_list, iado);

  ddos1[iado].setZero();

  int pos = g_index_list(0) - 1;
  if (pos != -1) {
    ddos1[iado] = -gams * total[pos];
    ddos1[iado] -= d->i * (total[pos] * d->hamt - d->hamt * total[pos]);
  }

  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    pos = g_index_list(i_nind + 1) - 1;
    if (pos != -1) {
      if (key(i_nind)) {
        ddos1[iado] -= (parity_key(i_nind) ? 1.0 : -1.0) * d->i *
                       (parity_p_n * total[pos] * d->qmdta_l[i_nind] -
                        d->qmdta_r[i_nind] * total[pos]);
      } else {
        ddos1[iado] -= (parity_key(i_nind) ? 1.0 : -1.0) * d->i *
                       (parity_p_n * total[pos] * d->qmdtc_l[i_nind] +
                        d->qmdtc_r[i_nind] * total[pos]);
      }
    }
  }

  init_nNNmat(ddos1[iado]);
}

template <typename FUN>
void syl_cal_ddos_gen(nNNmat &ddos1, const nNNmat &total,
                      const nNNmat &ddos2, const DEOM_DATA *d,
                      const DEOM_SYL *syl, const int iado,
                      FUN gen_g_index) {
  key_vec key(d->keys[iado]), parity_key(d->nind);
  double parity_p_n = 0;
  gen_parity(d, key, parity_key, &parity_p_n);
  ArrayXi g_index_list(d->nind + 1);
  gen_g_index(d, key, g_index_list, iado);

  ddos1[iado].setZero();

  int pos = g_index_list(0) - 1;
  if (pos != -1) ddos1[iado] = syl->OMG * total[pos] + ddos2[pos];

  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    pos = g_index_list(i_nind + 1) - 1;
    if ((pos != -1)) {
      if (key(i_nind)) {
        ddos1[iado] -= (parity_key(i_nind) ? 1.0 : -1.0) * d->i *
                       (parity_p_n * d->qmdtc_l[i_nind] * total[pos] +
                        total[pos] * d->qmdtc_r[i_nind]);
      } else {
        ddos1[iado] -= (parity_key(i_nind) ? 1.0 : -1.0) * d->i *
                       (parity_p_n * d->qmdta_l[i_nind] * total[pos] -
                        total[pos] * d->qmdta_r[i_nind]);
      }
    }
  }

  init_nNNmat(ddos1[iado]);

  ArrayXi key_syl(syl->nind);
  gen_syl_key(key_syl, key, syl);
  pos = generate_hash_value_syl(key_syl, d, syl) - 1;
  syl_complex_storge(ddos1[iado], syl, pos, d->ferr, ddos1[iado]);
}

template <typename FUN>
void syl_cal_gen(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
                 const DEOM_SYL *syl, const int iado, FUN gen_g_index) {
  key_vec key(d->keys[iado]), parity_key(d->nind);
  double parity_p_n = 0;
  gen_parity(d, key, parity_key, &parity_p_n);
  ArrayXi g_index_list(d->nind + 1);
  gen_g_index(d, key, g_index_list, iado);

  ddos1[iado].setZero();

  int pos = g_index_list(0) - 1;
  if (pos != -1) ddos1[iado] = syl->OMG * total[pos];

  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    pos = g_index_list(i_nind + 1) - 1;
    if ((pos != -1)) {
      if (key(i_nind)) {
        ddos1[iado] -= (parity_key(i_nind) ? 1.0 : -1.0) * d->i *
                       (parity_p_n * d->qmdtc_l[i_nind] * total[pos] +
                        total[pos] * d->qmdtc_r[i_nind]);
      } else {
        ddos1[iado] -= (parity_key(i_nind) ? 1.0 : -1.0) * d->i *
                       (parity_p_n * d->qmdta_l[i_nind] * total[pos] -
                        total[pos] * d->qmdta_r[i_nind]);
      }
    }
  }

  init_nNNmat(ddos1[iado]);

  ArrayXi key_syl(syl->nind);
  gen_syl_key(key_syl, key, syl);
  pos = generate_hash_value_syl(key_syl, d, syl) - 1;
  syl_complex_storge(ddos1[iado], syl, pos, d->ferr, ddos1[iado]);
}

template <typename FUN>
void syl_cal_Hei_gen(nNNmat &ddos1, const nNNmat &total,
                     const nNNmat &ddos2, const DEOM_DATA *d,
                     const DEOM_SYL *syl, const int iado,
                     FUN gen_g_index) {
  key_vec key(d->keys[iado]), parity_key(d->nind);
  double parity_p_n = 0;
  gen_parity(d, key, parity_key, &parity_p_n);
  ArrayXi g_index_list(d->nind + 1);
  gen_g_index(d, key, g_index_list, iado);

  ddos1[iado].setZero();

  int pos = g_index_list(0) - 1;
  if (pos != -1) ddos1[iado] = -syl->OMG * total[pos] - ddos2[pos];

  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    pos = g_index_list(i_nind + 1) - 1;
    if ((pos != -1)) {
      if (key(i_nind)) {
        ddos1[iado] += (parity_key(i_nind) ? 1.0 : -1.0) * d->i *
                       (parity_p_n * total[pos] * d->qmdta_l[i_nind] -
                        d->qmdta_r[i_nind] * total[pos]);
      } else {
        ddos1[iado] += (parity_key(i_nind) ? 1.0 : -1.0) * d->i *
                       (parity_p_n * total[pos] * d->qmdtc_l[i_nind] +
                        d->qmdtc_r[i_nind] * total[pos]);
      }
    }
  }

  init_nNNmat(ddos1[iado]);

  ArrayXi key_syl(syl->nind);
  gen_syl_key(key_syl, key, syl);
  pos = generate_hash_value_syl(key_syl, d, syl) - 1;
  syl_complex_Hei_storge(ddos1[iado], syl, pos, d->ferr, ddos1[iado]);
}

template <typename FUN>
void syl_cal_1_gen(nNNmat &ddos1, const nNNmat &total, const nNNmat &ddos2,
                   const DEOM_DATA *d, const DEOM_SYL *syl, const int iado,
                   FUN gen_g_index) {
  key_vec key(d->keys[iado]), parity_key(d->nind);
  ArrayXi key_syl(syl->nind);
  gen_syl_key(key_syl, key, syl);
  double parity_p_n = 0;
  gen_parity(d, key, parity_key, &parity_p_n);
  ArrayXi g_index_list(d->nind + 1);
  gen_g_index(d, key, g_index_list, iado);

  ddos1[iado].setZero();

  int pos = g_index_list(0) - 1;
  if (pos != -1) ddos1[iado] += syl->OMG * total[pos] + ddos2[pos];

  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    pos = g_index_list(i_nind + 1) - 1;
    if ((pos != -1)) {
      if (key(i_nind)) {
        ddos1[iado] -= (parity_key(i_nind) ? 1.0 : -1.0) * d->i *
                       (parity_p_n * d->qmdtc_l[i_nind] * total[pos] -
                        total[pos] * d->qmdtc_r[i_nind]);
      } else {
        ddos1[iado] -= (parity_key(i_nind) ? 1.0 : -1.0) * d->i *
                       (parity_p_n * d->qmdta_l[i_nind] * total[pos] +
                        total[pos] * d->qmdta_r[i_nind]);
      }
    }
  }

  init_nNNmat(ddos1[iado]);
  pos = generate_hash_value_syl(key_syl, d, syl) - 1;
  syl_complex_storge(ddos1[iado], syl, pos, d->ferr, ddos1[iado]);
}

template <typename FUN>
void syl_cal_1_Hei_gen(nNNmat &ddos1, const nNNmat &total,
                       const nNNmat &ddos2, const DEOM_DATA *d,
                       const DEOM_SYL *syl, const int iado,
                       FUN gen_g_index) {
  key_vec key(d->keys[iado]), parity_key(d->nind);
  ArrayXi key_syl(syl->nind);
  gen_syl_key(key_syl, key, syl);
  ArrayXi g_index_list(d->nind + 1);
  gen_g_index(d, key, g_index_list, iado);

  double parity_p_n = 0;
  gen_parity(d, key, parity_key, &parity_p_n);

  ddos1[iado].setZero();

  int pos = g_index_list(0) - 1;
  if (pos != -1) ddos1[iado] = -syl->OMG * total[pos] - ddos2[pos];

  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    pos = g_index_list(i_nind + 1) - 1;
    if ((pos != -1)) {
      if (key(i_nind)) {
        ddos1[iado] += (parity_key(i_nind) ? 1.0 : -1.0) * d->i *
                       (parity_p_n * total[pos] * d->qmdta_l[i_nind] +
                        d->qmdta_r[i_nind] * total[pos]);
      } else {
        ddos1[iado] += (parity_key(i_nind) ? 1.0 : -1.0) * d->i *
                       (parity_p_n * total[pos] * d->qmdtc_l[i_nind] -
                        d->qmdtc_r[i_nind] * total[pos]);
      }
    }
  }

  init_nNNmat(ddos1[iado]);
  pos = generate_hash_value_syl(key_syl, d, syl) - 1;
  syl_complex_Hei_storge(ddos1[iado], syl, pos, d->ferr, ddos1[iado]);
}

template <typename FUN>
void rem_cal_imag_gen(nNNmat &ddos1, const nNNmat &total,
                      const DEOM_DATA *d, const int iado,
                      FUN gen_g_index) {
  key_vec key(d->keys[iado]), parity_key(d->nind);
  complex<double> gams = 0;
  double parity_p_n = 0;
  gen_parity(d, key, parity_key, &parity_p_n, &gams);
  ArrayXi g_index_list(d->nind + 1);
  gen_g_index(d, key, g_index_list, iado);

  ddos1[iado].setZero();

  int pos = g_index_list(0) - 1;
  if (pos != -1) {
    ddos1[iado] = d->i * (gams * total[pos]);
    ddos1[iado] -= (d->hamt * total[pos] - total[pos] * d->hamt);
  }

  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    pos = g_index_list(i_nind + 1) - 1;
    if (pos != -1) {
      if (key(i_nind)) {
        ddos1[iado] -= (parity_key(i_nind) ? 1.0 : -1.0) *
                       (parity_p_n * d->qmdtc_l[i_nind] * total[pos]);
      } else {
        ddos1[iado] -= (parity_key(i_nind) ? 1.0 : -1.0) *
                       (parity_p_n * d->qmdta_l[i_nind] * total[pos]);
      }
    }
  }

  init_nNNmat(ddos1[iado]);
}

template <typename FUN>
void rem_cal_cfw_gen(nNNmat &ddos1, const nNNmat &total,
                     const DEOM_DATA *d, const int iado, FUN gen_g_index) {
  key_vec key(d->keys[iado]), parity_key(d->nind);
  complex<double> gams = 0;
  double parity_p_n = 0;
  gen_parity(d, key, parity_key, &parity_p_n, &gams);
  ArrayXi g_index_list(d->nind + 1);
  gen_g_index(d, key, g_index_list, iado);

  ddos1[iado].setZero();

  int pos = g_index_list(0) - 1;
  if (pos != -1) {
    ddos1[iado] = -(gams * total[pos]);
    ddos1[iado] -= d->i * (d->hamt * total[pos] - total[pos] * d->hamt);
  }

  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    pos = g_index_list(i_nind + 1) - 1;
    if (pos != -1) {
      if (key(i_nind)) {
        ddos1[iado] -= (parity_key(i_nind) ? 1.0 : -1.0) * d->i *
                       (d->cfw.Lambda_n * parity_p_n * d->qmdtc_l[i_nind] *
                            total[pos] +
                        total[pos] * d->cfw.Lambda_p * d->qmdtc_r[i_nind]);
      } else {
        ddos1[iado] -= (parity_key(i_nind) ? 1.0 : -1.0) * d->i *
                       (d->cfw.Lambda_n * parity_p_n * d->qmdta_l[i_nind] *
                            total[pos] -
                        total[pos] * d->cfw.Lambda_p * d->qmdta_r[i_nind]);
      }
    }
  }

  init_nNNmat(ddos1[iado]);
}