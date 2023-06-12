#pragma once
#include "algebra.hpp"
#include "deom.hpp"
#include "linearalgebra.hpp"

void construct_Mapping_filter_p(DEOM_DATA* d, Trie* tree, int iado) {
  key_vec key(d->keys[iado]), key_swell(d->nind);
  llint hash_val = gen_hash_value(key, d);
  int pos_get = tree->find(hash_val);
  d->g_index_list_pos[pos_get](0) = pos_get + 1;

  for (int mp = 0; mp < d->nind; mp++) {
    if (key.sum() < d->lmax) {
      hash_val = generate_key_plus(key, key_swell, mp, d);
      set_ddos(key_swell, mp + 1, pos_get, hash_val, tree, d);
    }

    if (key(mp) > 0) {
      hash_val = generate_key_minus(key, key_swell, mp, d);
      set_ddos(key_swell, mp + d->nind + 1, pos_get, hash_val, tree, d);
    }
  }
}

void construct_Mapping_filter_p(nNNmat& ddos, DEOM_DATA* d, Trie* tree,
                                int iado) {
  key_vec key(d->keys[iado]), key_swell(d->nind);
  llint hash_val = gen_hash_value(key, d);
  bool flag = is_valid(key, ddos[iado], d);

  if (flag || hash_val == 1) {
    int pos_get = tree->find(hash_val);
    d->g_index_list_pos[pos_get](0) = pos_get + 1;

    for (int mp = 0; mp < d->nind; mp++) {
      if (key.sum() < d->lmax) {
        hash_val = generate_key_plus(key, key_swell, mp, d);
        set_ddos(key_swell, mp + 1, pos_get, hash_val, tree, d);
      }

      if (key(mp) > 0) {
        hash_val = generate_key_minus(key, key_swell, mp, d);
        set_ddos(key_swell, mp + d->nind + 1, pos_get, hash_val, tree, d);
      }
    }
  }
}

void rem_oprt_filter_p(nNNmat& ddos1, const nNNmat& total,
                       const DEOM_DATA* d, const int iado,
                       const char lcr) {
  key_vec key(d->keys[iado]);
  int pos = d->g_index_list_pos[iado](0) - 1;

  ddos1[iado].setZero();

  if (pos != -1) {
    for (int mp = 0; mp < d->nmod; mp++) {
      if (lcr == 'l') {
        ddos1[iado] += d->dipole.sdip[mp] * total[pos];
      } else if (lcr == 'r') {
        ddos1[iado] += total[pos] * d->dipole.sdip[mp];
      } else if (lcr == 'c') {
        ddos1[iado] += d->dipole.sdip[mp] * total[pos] -
                       total[pos] * d->dipole.sdip[mp];
      }
    }
  }

  for (int mp = 0; mp < d->nind; mp++) {
    double n = key(mp);
    pos = d->g_index_list_pos[iado](mp + 1) - 1;
    if (pos != -1) {
      if (lcr == 'l') {
        ddos1[iado] += d->dipole.bdip1[mp] * sqrt(n) *
                       d->dipole.qmdtc_l[mp] * total[pos];
      } else if (lcr == 'r') {
        ddos1[iado] += d->dipole.bdip1[mp] * sqrt(n) * total[pos] *
                       d->dipole.qmdtc_r[mp];
      } else if (lcr == 'c') {
        ddos1[iado] += d->dipole.bdip1[mp] * sqrt(n) *
                       (d->dipole.qmdtc_l[mp] * total[pos] -
                        total[pos] * d->dipole.qmdtc_r[mp]);
      }
    }

    pos = d->g_index_list_pos[iado](mp + 1 + d->nind) - 1;
    if (pos != -1) {
      if (lcr == 'l') {
        ddos1[iado] += d->dipole.bdip1[mp] * sqrt(n + 1) *
                       d->dipole.qmdta_l[mp] * total[pos];
      } else if (lcr == 'r') {
        ddos1[iado] += d->dipole.bdip1[mp] * sqrt(n + 1) * total[pos] *
                       d->dipole.qmdta_r[mp];
      } else if (lcr == 'c') {
        ddos1[iado] += d->dipole.bdip1[mp] * sqrt(n + 1) *
                       (d->dipole.qmdta_l[mp] * total[pos] -
                        total[pos] * d->dipole.qmdta_r[mp]);
      }
    }
  }
}

void rem_cal_filter_p(nNNmat& ddos1, const nNNmat& total,
                      const DEOM_DATA* d, const int iado) {
  key_vec key(d->keys[iado]);
  complex<double> gams = 0.0;
  for (int i = 0; i < d->nind; i++) gams += (double)(key(i)) * d->expn(i);
  int pos = d->g_index_list_pos[iado](0) - 1;

  ddos1[iado].setZero();

  if (pos != -1) {
    ddos1[iado] = -gams * total[pos];
    ddos1[iado] -= d->i * (d->hamt * total[pos] - total[pos] * d->hamt);
  }

  for (int mp = 0; mp < d->nind; mp++) {
    double n = key(mp);
    pos = d->g_index_list_pos[iado](mp + 1) - 1;
    if (pos != -1)
      ddos1[iado] -=
          d->i * sqrt(n) *
          (d->qmdtc_l[mp] * total[pos] - total[pos] * d->qmdtc_r[mp]);

    pos = d->g_index_list_pos[iado](mp + 1 + d->nind) - 1;
    if (pos != -1)
      ddos1[iado] -=
          d->i * sqrt(n + 1) *
          (d->qmdta_l[mp] * total[pos] - total[pos] * d->qmdta_r[mp]);
  }
}

void rem_cal_hei_filter_p(nNNmat& ddos1, const nNNmat& total,
                          const DEOM_DATA* d, const int iado) {}

void rem_cal_1_filter_p(nNNmat& ddos1, const nNNmat& total,
                        const DEOM_DATA* d, const int iado) {}

void rem_cal_1_hei_filter_p(nNNmat& ddos1, const nNNmat& total,
                            const DEOM_DATA* d, const int iado) {}

void syl_cal_ddos_filter_p(nNNmat& ddos1, const nNNmat& total,
                           const nNNmat& ddos2, const DEOM_DATA* d,
                           const DEOM_SYL* syl, const int iado) {
  key_vec key(d->keys[iado]);
  complex<double> gams = 0.0;
  for (int i = 0; i < d->nind; i++) gams += (double)(key(i)) * d->expn(i);
  int pos = d->g_index_list_pos[iado](0) - 1;

  ddos1[iado].setZero();

  if (pos != -1) ddos1[iado] = syl->OMG * total[pos] + gams * ddos2[pos];

  for (int mp = 0; mp < d->nind; mp++) {
    double n = key(mp);
    pos = d->g_index_list_pos[iado](mp + 1) - 1;
    if (pos != -1)
      ddos1[iado] +=
          d->i * sqrt(n) *
          (d->qmdtc_l[mp] * total[pos] - total[pos] * d->qmdtc_r[mp]);

    pos = d->g_index_list_pos[iado](mp + 1 + d->nind) - 1;
    if (pos != -1)
      ddos1[iado] +=
          d->i * sqrt(n + 1) *
          (d->qmdta_l[mp] * total[pos] - total[pos] * d->qmdta_r[mp]);
  }

  MatrixNcd hlft = d->i * d->hamt + gams * MatrixNcd::Identity();
  MatrixNcd hrht = d->i * d->hamt - syl->OMG * MatrixNcd::Identity();
  syl_complex(ddos1[iado], hlft, hrht, ddos1[iado]);
}

void syl_cal_filter_p(nNNmat& ddos1, const nNNmat& total,
                      const DEOM_DATA* d, const DEOM_SYL* syl,
                      const int iado) {
  key_vec key(d->keys[iado]);
  complex<double> gams = 0.0;
  for (int i = 0; i < d->nind; i++) gams += (double)(key(i)) * d->expn(i);
  int pos = d->g_index_list_pos[iado](0) - 1;

  ddos1[iado].setZero();

  if (pos != -1) ddos1[iado] = syl->OMG * total[pos];

  for (int mp = 0; mp < d->nind; mp++) {
    double n = key(mp);
    pos = d->g_index_list_pos[iado](mp + 1) - 1;
    if (pos != -1)
      ddos1[iado] +=
          d->i * sqrt(n) *
          (d->qmdtc_l[mp] * total[pos] - total[pos] * d->qmdtc_r[mp]);

    pos = d->g_index_list_pos[iado](mp + 1 + d->nind) - 1;
    if (pos != -1)
      ddos1[iado] +=
          d->i * sqrt(n + 1) *
          (d->qmdta_l[mp] * total[pos] - total[pos] * d->qmdta_r[mp]);
  }

  MatrixNcd hlft = d->i * d->hamt + gams * MatrixNcd::Identity();
  MatrixNcd hrht = d->i * d->hamt - syl->OMG * MatrixNcd::Identity();
  syl_complex(ddos1[iado], hlft, hrht, ddos1[iado]);
}

void syl_cal_hei_filter_p(nNNmat& ddos1, const nNNmat& total,
                          const nNNmat& ddos2, const DEOM_DATA* d,
                          const DEOM_SYL* syl, const int iado) {}

void syl_cal_1_filter_p(nNNmat& ddos1, const nNNmat& total,
                        const nNNmat& ddos2, const DEOM_DATA* d,
                        const DEOM_SYL* syl, const int iado) {}

void syl_cal_1_hei_filter_p(nNNmat& ddos1, const nNNmat& total,
                            const nNNmat& ddos2, const DEOM_DATA* d,
                            const DEOM_SYL* syl, const int iado) {}

void rem_cal_imag_filter_p(nNNmat& ddos1, const nNNmat& total,
                           const DEOM_DATA* d, const int iado) {}

void rem_cal_cfw_filter_p(nNNmat& ddos1, const nNNmat& total,
                          const DEOM_DATA* d, const int iado) {}

complex<double> caculate_trace(DEOM_DATA* d, Trie* tree) {
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