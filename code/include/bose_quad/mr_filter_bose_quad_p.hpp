#pragma once

#include "algebra.hpp"
#include "cfw.hpp"
#include "deom.hpp"
#include "linearalgebra.hpp"

void construct_Mapping_filter_p(DEOM_DATA *d, Trie *tree, int iado) {
  key_vec key(d->nind), key_swell(d->nind), key_swell1(d->nind);
  key.noalias() = d->keys[iado];
  llint hash_val = gen_hash_value(key, d);
  int pos_get = tree->find(hash_val);
  d->g_index_list_pos[pos_get](0) = pos_get + 1;

  for (int mp1 = 0; mp1 < d->nind; mp1++) {
    if (key.sum() < d->lmax) {
      hash_val = generate_key_plus(key, key_swell, mp1, d);
      set_ddos(key_swell, mp1 + 1, pos_get, hash_val, tree, d);

      for (int mp2 = 0; mp2 < d->nind; mp2++) {
        if (key.sum() < d->lmax - 1) {
          hash_val = generate_key_plus(key_swell, key_swell1, mp2, d);
          set_ddos(key_swell1, 2 * d->nind + 1 + mp1 * d->nind + mp2,
                   pos_get, hash_val, tree, d);
        }

        if (key_swell(mp2) > 0) {
          hash_val = generate_key_minus(key_swell, key_swell1, mp2, d);
          set_ddos(
              key_swell1,
              2 * d->nind + 1 + d->nind * d->nind + mp1 * d->nind + mp2,
              pos_get, hash_val, tree, d);
        }
      }
    }

    if (key(mp1) > 0) {
      hash_val = generate_key_minus(key, key_swell, mp1, d);
      set_ddos(key_swell, mp1 + d->nind + 1, pos_get, hash_val, tree, d);

      for (int mp2 = 0; mp2 < d->nind; mp2++) {
        if (key_swell(mp2) > 0) {
          hash_val = generate_key_minus(key_swell, key_swell1, mp2, d);
          set_ddos(key_swell1,
                   2 * d->nind + 1 + 2 * d->nind * d->nind +
                       mp1 * d->nind + mp2,
                   pos_get, hash_val, tree, d);
        }
      }
    }
  }
}

void construct_Mapping_filter_p(nNNmat &ddos, DEOM_DATA *d, Trie *tree,
                                int iado) {
  key_vec key(d->keys[iado]), key_swell(d->nind), key_swell1(d->nind);
  int pos_get;
  llint hash_val = gen_hash_value(key, d);
  bool flag = is_valid(key, ddos[iado], d);
  if (flag || hash_val == 1) {
    pos_get = tree->find(hash_val);
    d->g_index_list_pos[pos_get](0) = pos_get + 1;

    for (int mp1 = 0; mp1 < d->nind; mp1++) {
      if (key.sum() < d->lmax) {
        hash_val = generate_key_plus(key, key_swell, mp1, d);
        set_ddos(key_swell, mp1 + 1, pos_get, hash_val, tree, d);

        for (int mp2 = 0; mp2 < d->nind; mp2++) {
          if (key.sum() < d->lmax - 1) {
            hash_val = generate_key_plus(key_swell, key_swell1, mp2, d);
            set_ddos(key_swell1, 2 * d->nind + 1 + mp1 * d->nind + mp2,
                     pos_get, hash_val, tree, d);
          }

          if (key_swell(mp2) > 0) {
            hash_val = generate_key_minus(key_swell, key_swell1, mp2, d);
            set_ddos(
                key_swell1,
                2 * d->nind + 1 + d->nind * d->nind + mp1 * d->nind + mp2,
                pos_get, hash_val, tree, d);
          }
        }
      }

      if (key(mp1) > 0) {
        hash_val = generate_key_minus(key, key_swell, mp1, d);
        set_ddos(key_swell, mp1 + d->nind + 1, pos_get, hash_val, tree, d);

        for (int mp2 = 0; mp2 < d->nind; mp2++) {
          if (key_swell(mp2) > 0) {
            hash_val = generate_key_minus(key_swell, key_swell1, mp2, d);
            set_ddos(key_swell1,
                     2 * d->nind + 1 + 2 * d->nind * d->nind +
                         mp1 * d->nind + mp2,
                     pos_get, hash_val, tree, d);
          }
        }
      }
    }
  }
}

void rem_oprt_filter_p(nNNmat &ddos1, const nNNmat &total,
                       const DEOM_DATA *d, const int iado,
                       const char lcr) {
  int pos = d->g_index_list_pos[iado](0) - 1;
  double n1, n12;
  key_vec key(d->keys[iado]);
  ddos1[iado].setZero();

  if (pos != -1) {
    for (int nmod = 0; nmod < d->nmod; nmod++) {
      if (lcr == 'l') {
        ddos1[iado].noalias() +=
            (d->dipole.sdip[nmod] + d->dipole.renormalize[nmod]) *
            total[pos];
      } else if (lcr == 'r') {
        ddos1[iado].noalias() +=
            total[pos] *
            (d->dipole.sdip[nmod] + d->dipole.renormalize[nmod]);
      } else if (lcr == 'c') {
        ddos1[iado].noalias() +=
            (d->dipole.sdip[nmod] + d->dipole.renormalize[nmod]) *
                total[pos] -
            total[pos] *
                (d->dipole.sdip[nmod] + d->dipole.renormalize[nmod]);
      }
    }
  }

  for (int mp1 = 0; mp1 < d->nind; mp1++) {
    pos = d->g_index_list_pos[iado](mp1 + 1) - 1;
    n1 = key(mp1);
    if (pos != -1) {
      if (lcr == 'l') {
        ddos1[iado].noalias() += d->dipole.bdip1[mp1] * d->alp1 *
                                 sqrt(n1) * d->dipole.qmdtc_l[mp1] *
                                 total[pos];
      } else if (lcr == 'r') {
        ddos1[iado].noalias() += d->dipole.bdip1[mp1] * d->alp1 *
                                 sqrt(n1) * total[pos] *
                                 d->dipole.qmdtc_r[mp1];
      } else if (lcr == 'c') {
        ddos1[iado].noalias() += d->dipole.bdip1[mp1] * d->alp1 *
                                 sqrt(n1) *
                                 (d->dipole.qmdtc_l[mp1] * total[pos] -
                                  total[pos] * d->dipole.qmdtc_r[mp1]);
      }
    }

    pos = d->g_index_list_pos[iado](mp1 + 1 + d->nind) - 1;
    if (pos != -1) {
      if (lcr == 'l') {
        ddos1[iado].noalias() += d->dipole.bdip1[mp1] * d->alp1 *
                                 sqrt(n1 + 1) * d->dipole.qmdta_l[mp1] *
                                 total[pos];
      } else if (lcr == 'r') {
        ddos1[iado].noalias() += d->dipole.bdip1[mp1] * d->alp1 *
                                 sqrt(n1 + 1) * total[pos] *
                                 d->dipole.qmdta_r[mp1];
      } else if (lcr == 'c') {
        ddos1[iado].noalias() += d->dipole.bdip1[mp1] * d->alp1 *
                                 sqrt(n1 + 1) *
                                 (d->dipole.qmdta_l[mp1] * total[pos] -
                                  total[pos] * d->dipole.qmdta_r[mp1]);
      }
    }

    for (int mp2 = 0; mp2 < d->nind; mp2++) {
      pos = d->g_index_list_pos[iado](2 * d->nind + 1 + mp1 * d->nind +
                                      mp2) -
            1;
      if (pos != -1) {
        n12 = (mp1 == mp2) ? (key(mp1) * (key(mp2) - 1.0))
                           : (key(mp1) * key(mp2));
        if (lcr == 'l') {
          ddos1[iado].noalias() +=
              d->dipole.bdip2[mp1] * d->alp2 * sqrt(n12) *
              d->dipole.qmdt2c_l[INDEX2(mp1, mp2, d->nind)] * total[pos];
        } else if (lcr == 'r') {
          ddos1[iado].noalias() +=
              d->dipole.bdip2[mp1] * d->alp2 * sqrt(n12) * total[pos] *
              d->dipole.qmdt2c_r[INDEX2(mp1, mp2, d->nind)];
        } else if (lcr == 'c') {
          ddos1[iado].noalias() +=
              d->dipole.bdip2[mp1] * d->alp2 * sqrt(n12) *
              (d->dipole.qmdt2c_l[INDEX2(mp1, mp2, d->nind)] * total[pos] -
               total[pos] * d->dipole.qmdt2c_r[INDEX2(mp1, mp2, d->nind)]);
        }
      }

      pos = d->g_index_list_pos[iado](2 * d->nind + 1 + d->nind * d->nind +
                                      mp1 * d->nind + mp2) -
            1;
      if (pos != -1) {
        n12 = (mp1 == mp2) ? (key(mp1) * key(mp2))
                           : (key(mp1) * (key(mp2) + 1.0));
        if (lcr == 'l') {
          ddos1[iado].noalias() +=
              2.0 * d->dipole.bdip2[mp1] * d->alp2 * sqrt(n12) *
              d->dipole.qmdt2b_l[INDEX2(mp1, mp2, d->nind)] * total[pos];
        } else if (lcr == 'r') {
          ddos1[iado].noalias() +=
              2.0 * d->dipole.bdip2[mp1] * d->alp2 * sqrt(n12) *
              total[pos] * d->dipole.qmdt2b_r[INDEX2(mp1, mp2, d->nind)];
        } else if (lcr == 'c') {
          ddos1[iado].noalias() +=
              2.0 * d->dipole.bdip2[mp1] * d->alp2 * sqrt(n12) *
              (d->dipole.qmdt2b_l[INDEX2(mp1, mp2, d->nind)] * total[pos] -
               total[pos] * d->dipole.qmdt2b_r[INDEX2(mp1, mp2, d->nind)]);
        }
      }

      pos = d->g_index_list_pos[iado](2 * d->nind + 1 +
                                      2 * d->nind * d->nind +
                                      mp1 * d->nind + mp2) -
            1;
      if (pos != -1) {
        n12 = (mp1 == mp2) ? ((key(mp1) + 1.0) * (key(mp2) + 2.0))
                           : ((key(mp1) + 1.0) * (key(mp2) + 1.0));
        if (lcr == 'l') {
          ddos1[iado].noalias() +=
              d->dipole.bdip2[mp1] * d->alp2 * sqrt(n12) *
              d->dipole.qmdt2a_l[INDEX2(mp1, mp2, d->nind)] * total[pos];
        } else if (lcr == 'r') {
          ddos1[iado].noalias() +=
              d->dipole.bdip2[mp1] * d->alp2 * sqrt(n12) * total[pos] *
              d->dipole.qmdt2a_r[INDEX2(mp1, mp2, d->nind)];
        } else if (lcr == 'c') {
          ddos1[iado].noalias() +=
              d->dipole.bdip2[mp1] * d->alp2 * sqrt(n12) *
              (d->dipole.qmdt2a_l[INDEX2(mp1, mp2, d->nind)] * total[pos] -
               total[pos] * d->dipole.qmdt2a_r[INDEX2(mp1, mp2, d->nind)]);
        }
      }
    }
  }
}

void rem_oprt_hei_filter_p(nNNmat &ddos1, const nNNmat &total,
                           const DEOM_DATA *d, const int iado,
                           const char lcr) {}

void rem_cal_filter_p(nNNmat &ddos1, const nNNmat &total,
                      const DEOM_DATA *d, const int iado) {
  ddos1[iado].setZero();
  int pos = d->g_index_list_pos[iado](0) - 1;
  key_vec key(d->keys[iado]);
  complex<double> gams = 0.0;
  for (int i = 0; i < d->nind; i++) gams += (double)(key(i)) * d->expn(i);

  if (pos != -1) {
    ddos1[iado].noalias() -= gams * total[pos];
    ddos1[iado].noalias() -=
        d->i * (d->hamt * total[pos] - total[pos] * d->hamt);
    for (int nmod = 0; nmod < d->nmod; nmod++)
      ddos1[iado].noalias() -= d->i * (d->renormalize[nmod] * total[pos] -
                                       total[pos] * d->renormalize[nmod]);
  }

  for (int mp1 = 0; mp1 < d->nind; mp1++) {
    pos = d->g_index_list_pos[iado](mp1 + 1) - 1;
    double n1 = key(mp1);
    if (pos != -1) {
      ddos1[iado].noalias() -=
          d->alp1 * d->i * sqrt(n1) *
          (d->qmdtc_l[mp1] * total[pos] - total[pos] * d->qmdtc_r[mp1]);
    }

    pos = d->g_index_list_pos[iado](mp1 + 1 + d->nind) - 1;
    if (pos != -1) {
      ddos1[iado].noalias() -=
          d->alp1 * d->i * sqrt(n1 + 1) *
          (d->qmdta_l[mp1] * total[pos] - total[pos] * d->qmdta_r[mp1]);
    }

    for (int mp2 = 0; mp2 < d->nind; mp2++) {
      pos = d->g_index_list_pos[iado](2 * d->nind + 1 + mp1 * d->nind +
                                      mp2) -
            1;
      if (pos != -1) {
        double n12 = (mp1 == mp2) ? (key(mp1) * (key(mp2) - 1.0))
                                  : (key(mp1) * key(mp2));
        ddos1[iado].noalias() -=
            d->alp2 * d->i * sqrt(n12) *
            (d->qmdt2c_l[INDEX2(mp1, mp2, d->nind)] * total[pos] -
             total[pos] * d->qmdt2c_r[INDEX2(mp1, mp2, d->nind)]);
      }

      pos = d->g_index_list_pos[iado](2 * d->nind + 1 + d->nind * d->nind +
                                      mp1 * d->nind + mp2) -
            1;
      if (pos != -1) {
        double n12 = (mp1 == mp2) ? (key(mp1) * key(mp2))
                                  : (key(mp1) * (key(mp2) + 1.0));
        ddos1[iado].noalias() -=
            2 * d->alp2 * d->i * sqrt(n12) *
            (d->qmdt2b_l[INDEX2(mp1, mp2, d->nind)] * total[pos] -
             total[pos] * d->qmdt2b_r[INDEX2(mp1, mp2, d->nind)]);
      }

      pos = d->g_index_list_pos[iado](2 * d->nind + 1 +
                                      2 * d->nind * d->nind +
                                      mp1 * d->nind + mp2) -
            1;
      if (pos != -1) {
        double n12 = (mp1 == mp2) ? ((key(mp1) + 1.0) * (key(mp2) + 2.0))
                                  : ((key(mp1) + 1.0) * (key(mp2) + 1.0));
        ddos1[iado].noalias() -=
            d->alp2 * d->i * sqrt(n12) *
            (d->qmdt2a_l[INDEX2(mp1, mp2, d->nind)] * total[pos] -
             total[pos] * d->qmdt2a_r[INDEX2(mp1, mp2, d->nind)]);
      }
    }
  }
}

void rem_cal_hei_filter_p(nNNmat &ddos1, const nNNmat &total,
                          const DEOM_DATA *d, const int iado) {}

void rem_cal_1_filter_p(nNNmat &ddos1, const nNNmat &total,
                        const DEOM_DATA *d, const int iado) {}

void rem_cal_1_hei_filter_p(nNNmat &ddos1, const nNNmat &total,
                            const DEOM_DATA *d, const int iado) {}

void syl_cal_ddos_filter_p(nNNmat &ddos1, const nNNmat &total,
                           const nNNmat &ddos2, const DEOM_DATA *d,
                           const DEOM_SYL *syl, const int iado) {}

void syl_cal_filter_p(nNNmat &ddos1, const nNNmat &total,
                      const DEOM_DATA *d, const DEOM_SYL *syl,
                      const int iado) {}

void syl_cal_hei_filter_p(nNNmat &ddos1, const nNNmat &total,
                          const nNNmat &ddos2, const DEOM_DATA *d,
                          const DEOM_SYL *syl, const int iado) {}

void syl_cal_1_filter_p(nNNmat &ddos1, const nNNmat &total,
                        const nNNmat &ddos2, const DEOM_DATA *d,
                        const DEOM_SYL *syl, const int iado) {}

void syl_cal_1_hei_filter_p(nNNmat &ddos1, const nNNmat &total,
                            const nNNmat &ddos2, const DEOM_DATA *d,
                            const DEOM_SYL *syl, const int iado) {}

void rem_cal_imag_filter_p(nNNmat &ddos1, const nNNmat &total,
                           const DEOM_DATA *d, const int iado) {
  ddos1[iado].setZero();
  int pos = d->g_index_list_pos[iado](0) - 1;
  key_vec key(d->keys[iado]);
  complex<double> gams = 0.0;
  for (int i = 0; i < d->nind; i++) gams += (double)(key(i)) * d->expn(i);

  if (pos != -1) {
    ddos1[iado] -= -d->i * gams * total[pos];
    ddos1[iado] -= d->hamt * total[pos] - total[pos] * d->hamt;
    for (int nmod = 0; nmod < d->nmod; nmod++)
      ddos1[iado].noalias() -= d->renormalize[nmod] * total[pos];
  }

  for (int mp1 = 0; mp1 < d->nind; mp1++) {
    pos = d->g_index_list_pos[iado](mp1 + 1) - 1;
    double n1 = key(mp1);
    if (pos != -1) {
      ddos1[iado] -= d->alp1 * sqrt(n1) * d->qmdtc_l[mp1] * total[pos];
    }

    pos = d->g_index_list_pos[iado](mp1 + 1 + d->nind) - 1;
    if (pos != -1) {
      ddos1[iado] -= d->alp1 * sqrt(n1 + 1) * d->qmdta_l[mp1] * total[pos];
    }

    for (int mp2 = 0; mp2 < d->nind; mp2++) {
      pos = d->g_index_list_pos[iado](2 * d->nind + 1 + mp1 * d->nind +
                                      mp2) -
            1;
      if (pos != -1) {
        double n12 = (mp1 == mp2) ? (key(mp1) * (key(mp2) - 1.0))
                                  : (key(mp1) * key(mp2));
        ddos1[iado] -= d->alp2 * sqrt(n12) *
                       d->qmdt2c_l[INDEX2(mp1, mp2, d->nind)] * total[pos];
      }

      pos = d->g_index_list_pos[iado](2 * d->nind + 1 + d->nind * d->nind +
                                      mp1 * d->nind + mp2) -
            1;
      if (pos != -1) {
        double n12 = (mp1 == mp2) ? (key(mp1) * key(mp2))
                                  : (key(mp1) * (key(mp2) + 1.0));
        ddos1[iado] -= 2 * d->alp2 * sqrt(n12) *
                       d->qmdt2b_l[INDEX2(mp1, mp2, d->nind)] * total[pos];
      }

      pos = d->g_index_list_pos[iado](2 * d->nind + 1 +
                                      2 * d->nind * d->nind +
                                      mp1 * d->nind + mp2) -
            1;
      if (pos != -1) {
        double n12 = (mp1 == mp2) ? ((key(mp1) + 1.0) * (key(mp2) + 2.0))
                                  : ((key(mp1) + 1.0) * (key(mp2) + 1.0));
        ddos1[iado] -= d->alp2 * sqrt(n12) *
                       d->qmdt2a_l[INDEX2(mp1, mp2, d->nind)] * total[pos];
      }
    }
  }
}

void rem_cal_cfw_filter_p(nNNmat &ddos1, const nNNmat &total,
                          const DEOM_DATA *d, const int iado) {
  int pos = d->g_index_list_pos[iado](0) - 1;
  ddos1[iado].setZero();
  key_vec key(d->keys[iado]);
  complex<double> gams = 0.0;
  for (int i = 0; i < d->nind; i++) gams += (double)(key(i)) * d->expn(i);

  if (pos != -1) {
    ddos1[iado] -= gams * total[pos];
    ddos1[iado].noalias() -=
        d->i * (d->hamt * total[pos] - total[pos] * d->hamt);
    for (int nmod = 0; nmod < d->nmod; nmod++) {
      ddos1[iado].noalias() -=
          d->i * (d->cfw.Lambda_n * d->renormalize[nmod] * total[pos] -
                  total[pos] * d->cfw.Lambda_p * d->renormalize[nmod]);
    }
  }

  for (int mp1 = 0; mp1 < d->nind; mp1++) {
    pos = d->g_index_list_pos[iado](mp1 + 1) - 1;
    double n1 = key(mp1);
    if (pos != -1) {
      ddos1[iado] -= d->alp1 * d->i * sqrt(n1) *
                     (d->cfw.Lambda_n * d->qmdtc_l[mp1] * total[pos] -
                      total[pos] * d->cfw.Lambda_p * d->qmdtc_r[mp1]);
    }

    pos = d->g_index_list_pos[iado](mp1 + 1 + d->nind) - 1;
    if (pos != -1) {
      ddos1[iado] -= d->alp1 * d->i * sqrt(n1 + 1) *
                     (d->cfw.Lambda_n * d->qmdta_l[mp1] * total[pos] -
                      total[pos] * d->cfw.Lambda_p * d->qmdta_r[mp1]);
    }

    for (int mp2 = 0; mp2 < d->nind; mp2++) {
      pos = d->g_index_list_pos[iado](2 * d->nind + 1 + mp1 * d->nind +
                                      mp2) -
            1;
      if (pos != -1) {
        double n12 = (mp1 == mp2) ? (key(mp1) * (key(mp2) - 1.0))
                                  : (key(mp1) * key(mp2));
        ddos1[iado] -=
            d->alp2 * d->i * sqrt(n12) *
            (d->cfw.Lambda_n * d->qmdt2c_l[INDEX2(mp1, mp2, d->nind)] *
                 total[pos] -
             total[pos] * d->cfw.Lambda_p *
                 d->qmdt2c_r[INDEX2(mp1, mp2, d->nind)]);
      }

      pos = d->g_index_list_pos[iado](2 * d->nind + 1 + d->nind * d->nind +
                                      mp1 * d->nind + mp2) -
            1;
      if (pos != -1) {
        double n12 = (mp1 == mp2) ? (key(mp1) * key(mp2))
                                  : (key(mp1) * (key(mp2) + 1.0));
        ddos1[iado] -=
            2 * d->alp2 * d->i * sqrt(n12) *
            (d->cfw.Lambda_n * d->qmdt2b_l[INDEX2(mp1, mp2, d->nind)] *
                 total[pos] -
             total[pos] * d->cfw.Lambda_p *
                 d->qmdt2b_r[INDEX2(mp1, mp2, d->nind)]);
      }

      pos = d->g_index_list_pos[iado](2 * d->nind + 1 +
                                      2 * d->nind * d->nind +
                                      mp1 * d->nind + mp2) -
            1;
      if (pos != -1) {
        double n12 = (mp1 == mp2) ? ((key(mp1) + 1.0) * (key(mp2) + 2.0))
                                  : ((key(mp1) + 1.0) * (key(mp2) + 1.0));
        ddos1[iado] -=
            d->alp2 * d->i * sqrt(n12) *
            (d->cfw.Lambda_n * d->qmdt2a_l[INDEX2(mp1, mp2, d->nind)] *
                 total[pos] -
             total[pos] * d->cfw.Lambda_p *
                 d->qmdt2a_r[INDEX2(mp1, mp2, d->nind)]);
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

  if (pos != -1)
    for (int nmod = 0; nmod < d->nmod; nmod++)
      sum_Matrix += d->dipole1.renormalize[nmod] * d->ddos[pos];

  for (int mp1 = 0; mp1 < d->nind; mp1++) {
    llint hash_val = generate_key_plus(d->zerokey, d->emptykey, mp1, d);
    pos = tree->find(hash_val);
    if (pos != -1)
      sum_Matrix += d->dipole1.bdip1[mp1] * d->alp1 *
                    d->dipole1.qmdta_l[mp1] * d->ddos[pos];

    for (int mp2 = 0; mp2 < d->nind; mp2++) {
      llint hash_val = generate_hash_value_plus(d->emptykey, mp2, d);
      pos = tree->find(hash_val);
      if (pos != -1) {
        sum_Matrix += d->dipole1.bdip2[mp1] * d->alp2 *
                      sqrt((mp1 == mp2) ? 2.0 : 1.0) *
                      d->dipole1.qmdt2a_l[INDEX2(mp1, mp2, d->nind)] *
                      d->ddos[pos];
      }
    }
  }

  return sum_Matrix.trace();
}