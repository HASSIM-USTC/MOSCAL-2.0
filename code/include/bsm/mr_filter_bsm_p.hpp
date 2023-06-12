#pragma once
#include "algebra.hpp"
#include "deom.hpp"
#include "linearalgebra.hpp"

void construct_Mapping_filter_p(DEOM_DATA* d, Trie* tree, int iado) {
  key_vec key(d->keys[iado]), key_swell1(d->nind), key_swell2(d->nind);
  llint hash_val = gen_hash_value(key, d);
  int pos_get = tree->find(hash_val);
  d->g_index_list_pos[pos_get](0) = pos_get + 1;

  int key_sum = key.sum();
  int key_sum_fp = key(0) + key(1);
  int key_sum_ma = key_sum - key_sum_fp;

  if (key[1] > 0) {
    generate_key_plus(key, key_swell1, 0, d);
    hash_val = generate_key_minus(key_swell1, key_swell2, 1, d);
    set_ddos(key_swell2, 1, pos_get, hash_val, tree, d);
  }

  if (key[0] > 0) {
    generate_key_plus(key, key_swell1, 1, d);
    hash_val = generate_key_minus(key_swell1, key_swell2, 0, d);
    set_ddos(key_swell2, 2, pos_get, hash_val, tree, d);
  }

  if ((key_sum_fp < d->lmax_fp - 1) && (key_sum < d->lmax - 1)) {
    for (int mp = 0; mp < 2; ++mp) {
      generate_key_plus(key, key_swell1, mp, d);
      hash_val = generate_key_plus(key_swell1, key_swell2, mp, d);
      set_ddos(key_swell2, mp + 3, pos_get, hash_val, tree, d);
    }

    generate_key_plus(key, key_swell1, 0, d);
    hash_val = generate_key_plus(key_swell1, key_swell2, 1, d);
    set_ddos(key_swell2, 5, pos_get, hash_val, tree, d);
  }

  for (int mp = 0; mp < 2; mp++) {
    if (key(mp) > 0) {
      hash_val = generate_key_minus(key, key_swell2, mp, d);
      set_ddos(key_swell2, 6 + mp, pos_get, hash_val, tree, d);
    }
  }

  if ((key_sum_fp < d->lmax_fp) && (key_sum < d->lmax)) {
    for (int mp = 0; mp < 2; ++mp) {
      hash_val = generate_key_plus(key, key_swell2, mp, d);
      set_ddos(key_swell2, 8 + mp, pos_get, hash_val, tree, d);
    }
  }

  for (int mp = 0; mp < 2; ++mp) {
    if ((key(mp) > 1)) {
      generate_key_minus(key, key_swell1, mp, d);
      hash_val = generate_key_minus(key_swell1, key_swell2, mp, d);
      set_ddos(key_swell2, 10 + mp, pos_get, hash_val, tree, d);
    }
  }

  if ((key(0) > 0) && (key(1) > 0)) {
    generate_key_minus(key, key_swell1, 0, d);
    hash_val = generate_key_minus(key_swell1, key_swell2, 1, d);
    set_ddos(key_swell2, 12, pos_get, hash_val, tree, d);
  }

  for (int mp1 = 0; mp1 < 2; ++mp1) {
    for (int mp2 = 2; mp2 < d->nind; ++mp2) {
      if ((key(mp2) > 0) &&
          ((key_sum_fp < d->lmax_fp) && (key_sum < d->lmax))) {
        generate_key_plus(key, key_swell1, mp1, d);
        hash_val = generate_key_minus(key_swell1, key_swell2, mp2, d);
        int allocate_pos = 11 + (d->nind - 2) * mp1 + mp2;
        set_ddos(key_swell2, allocate_pos, pos_get, hash_val, tree, d);
      }
    }
  }

  for (int mp1 = 0; mp1 < 2; ++mp1) {
    for (int mp2 = 2; mp2 < d->nind; ++mp2) {
      if ((key(mp1) > 0) &&
          ((key_sum_ma < d->lmax_ma) && (key_sum < d->lmax))) {
        generate_key_plus(key, key_swell1, mp2, d);
        hash_val = generate_key_minus(key_swell1, key_swell2, mp1, d);
        int allocate_pos = 11 + (d->nind - 2) * (mp1 + 2) + mp2;
        set_ddos(key_swell2, allocate_pos, pos_get, hash_val, tree, d);
      }
    }
  }

  if ((key_sum_fp < d->lmax_fp) && (key_sum_ma < d->lmax_ma) &&
      (key_sum < d->lmax - 1)) {
    for (int mp1 = 0; mp1 < 2; ++mp1) {
      for (int mp2 = 2; mp2 < d->nind; ++mp2) {
        generate_key_plus(key, key_swell1, mp1, d);
        hash_val = generate_key_plus(key_swell1, key_swell2, mp2, d);
        int allocate_pos = 11 + (4 + mp1) * (d->nind - 2) + mp2;
        set_ddos(key_swell2, allocate_pos, pos_get, hash_val, tree, d);
      }
    }
  }
}

void construct_Mapping_filter_p(nNNmat& ddos, DEOM_DATA* d, Trie* tree,
                                int iado) {
  key_vec key(d->keys[iado]), key_swell1(d->nind), key_swell2(d->nind);
  llint hash_val = gen_hash_value(key, d);

  bool flag = is_valid(key, ddos[iado], d);

  if (flag || hash_val == 1) {
    int pos_get = tree->find(hash_val);
    d->g_index_list_pos[pos_get](0) = pos_get + 1;

    int key_sum = key.sum();
    int key_sum_fp = key(0) + key(1);
    int key_sum_ma = key_sum - key_sum_fp;

    if (key[1] > 0) {
      generate_key_plus(key, key_swell1, 0, d);
      hash_val = generate_key_minus(key_swell1, key_swell2, 1, d);
      set_ddos(key_swell2, 1, pos_get, hash_val, tree, d);
    }

    if (key[0] > 0) {
      generate_key_plus(key, key_swell1, 1, d);
      hash_val = generate_key_minus(key_swell1, key_swell2, 0, d);
      set_ddos(key_swell2, 2, pos_get, hash_val, tree, d);
    }

    if ((key_sum_fp < d->lmax_fp - 1) && (key_sum < d->lmax - 1)) {
      for (int mp = 0; mp < 2; ++mp) {
        generate_key_plus(key, key_swell1, mp, d);
        hash_val = generate_key_plus(key_swell1, key_swell2, mp, d);
        set_ddos(key_swell2, mp + 3, pos_get, hash_val, tree, d);
      }

      generate_key_plus(key, key_swell1, 0, d);
      hash_val = generate_key_plus(key_swell1, key_swell2, 1, d);
      set_ddos(key_swell2, 5, pos_get, hash_val, tree, d);
    }

    for (int mp = 0; mp < 2; mp++) {
      if (key(mp) > 0) {
        hash_val = generate_key_minus(key, key_swell2, mp, d);
        set_ddos(key_swell2, 6 + mp, pos_get, hash_val, tree, d);
      }
    }

    if ((key_sum_fp < d->lmax_fp) && (key_sum < d->lmax)) {
      for (int mp = 0; mp < 2; ++mp) {
        hash_val = generate_key_plus(key, key_swell2, mp, d);
        set_ddos(key_swell2, 8 + mp, pos_get, hash_val, tree, d);
      }
    }

    for (int mp = 0; mp < 2; ++mp) {
      if ((key(mp) > 1)) {
        generate_key_minus(key, key_swell1, mp, d);
        hash_val = generate_key_minus(key_swell1, key_swell2, mp, d);
        set_ddos(key_swell2, 10 + mp, pos_get, hash_val, tree, d);
      }
    }

    if ((key(0) > 0) && (key(1) > 0)) {
      generate_key_minus(key, key_swell1, 0, d);
      hash_val = generate_key_minus(key_swell1, key_swell2, 1, d);
      set_ddos(key_swell2, 12, pos_get, hash_val, tree, d);
    }

    for (int mp1 = 0; mp1 < 2; ++mp1) {
      for (int mp2 = 2; mp2 < d->nind; ++mp2) {
        if ((key(mp2) > 0) &&
            ((key_sum_fp < d->lmax_fp) && (key_sum < d->lmax))) {
          generate_key_plus(key, key_swell1, mp1, d);
          hash_val = generate_key_minus(key_swell1, key_swell2, mp2, d);
          int allocate_pos = 11 + (d->nind - 2) * mp1 + mp2;
          set_ddos(key_swell2, allocate_pos, pos_get, hash_val, tree, d);
        }
      }
    }

    for (int mp1 = 0; mp1 < 2; ++mp1) {
      for (int mp2 = 2; mp2 < d->nind; ++mp2) {
        if ((key(mp1) > 0) &&
            ((key_sum_ma < d->lmax_ma) && (key_sum < d->lmax))) {
          generate_key_plus(key, key_swell1, mp2, d);
          hash_val = generate_key_minus(key_swell1, key_swell2, mp1, d);
          int allocate_pos = 11 + (d->nind - 2) * (mp1 + 2) + mp2;
          set_ddos(key_swell2, allocate_pos, pos_get, hash_val, tree, d);
        }
      }
    }

    if ((key_sum_fp < d->lmax_fp) && (key_sum_ma < d->lmax_ma) &&
        (key_sum < d->lmax - 1)) {
      for (int mp1 = 0; mp1 < 2; ++mp1) {
        for (int mp2 = 2; mp2 < d->nind; ++mp2) {
          generate_key_plus(key, key_swell1, mp1, d);
          hash_val = generate_key_plus(key_swell1, key_swell2, mp2, d);
          int allocate_pos = 11 + (4 + mp1) * (d->nind - 2) + mp2;
          set_ddos(key_swell2, allocate_pos, pos_get, hash_val, tree, d);
        }
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

  if (pos != -1) {
    double n1 = key(0);
    double n2 = key(1);
    complex<double> cl =
        d->dipole.bdip2[0] * d->alp2 *
        ((2 * n1 + 1) * d->coef_lft[0] + (2 * n2 + 1) * d->coef_lft[1]);
    complex<double> cr =
        d->dipole.bdip2[0] * d->alp2 *
        ((2 * n1 + 1) * d->coef_rht[0] + (2 * n2 + 1) * d->coef_rht[1]);
    if (lcr == 'l') {
      ddos1[iado] +=
          d->alp0 * d->dipole.bdip0[0] * (d->dipole.pdip0 * total[pos]) +
          cl * d->dipole.bdip2[0] * d->dipole.pdip2 * total[pos];
    } else if (lcr == 'r') {
      ddos1[iado] +=
          d->alp0 * d->dipole.bdip0[0] * (total[pos] * d->dipole.pdip0) +
          cr * d->dipole.bdip2[0] * total[pos] * d->dipole.pdip2;
    } else if (lcr == 'c') {
      ddos1[iado] +=
          d->alp0 * d->dipole.bdip0[0] *
              (d->dipole.pdip0 * total[pos] -
               total[pos] * d->dipole.pdip0) +
          cl * d->dipole.bdip2[0] * d->dipole.pdip2 * total[pos] -
          cr * d->dipole.bdip2[0] * total[pos] * d->dipole.pdip2;
    }
  }

  pos = d->g_index_list_pos[iado](1) - 1;
  if (pos != -1) {
    double n1 = key(0);
    double n2 = key(1);
    complex<double> sn =
        sqrt(n1 * (n2 + 1)) / d->coef_abs[0] * d->coef_abs[1];
    complex<double> cl =
        d->dipole.bdip2[0] * d->alp2 * 2. * sn * d->coef_lft[0];
    complex<double> cr =
        d->dipole.bdip2[0] * d->alp2 * 2. * sn * d->coef_rht[0];
    if (lcr == 'l') {
      ddos1[iado] += cl * d->dipole.pdip2 * total[pos];
    } else if (lcr == 'r') {
      ddos1[iado] += cr * total[pos] * d->dipole.pdip2;
    } else if (lcr == 'c') {
      ddos1[iado] += cl * d->dipole.pdip2 * total[pos] -
                     cr * total[pos] * d->dipole.pdip2;
    }
  }

  pos = d->g_index_list_pos[iado](2) - 1;
  if (pos != -1) {
    double n1 = key(0);
    double n2 = key(1);
    complex<double> sn =
        sqrt((n1 + 1) * n2) * d->coef_abs[0] / d->coef_abs[1];
    complex<double> cl =
        d->dipole.bdip2[0] * d->alp2 * 2. * sn * d->coef_lft[1];
    complex<double> cr =
        d->dipole.bdip2[0] * d->alp2 * 2. * sn * d->coef_rht[1];
    if (lcr == 'l') {
      ddos1[iado] += cl * d->dipole.pdip2 * total[pos];
    } else if (lcr == 'r') {
      ddos1[iado] += cr * total[pos] * d->dipole.pdip2;
    } else if (lcr == 'c') {
      ddos1[iado] += cl * d->dipole.pdip2 * total[pos] -
                     cr * total[pos] * d->dipole.pdip2;
    }
  }

  for (int mp = 0; mp < 2; ++mp) {
    pos = d->g_index_list_pos[iado](mp + 3) - 1;
    if (pos != -1) {
      double n1 = key(mp);
      complex<double> sn =
          sqrt(n1 * (n1 - 1)) / d->coef_abs[mp] / d->coef_abs[mp];
      complex<double> cl = d->dipole.bdip2[0] * d->alp2 * sn *
                           d->coef_lft[mp] * d->coef_lft[mp];
      complex<double> cr = d->dipole.bdip2[0] * d->alp2 * sn *
                           d->coef_rht[mp] * d->coef_rht[mp];
      if (lcr == 'l') {
        ddos1[iado] += cl * d->dipole.pdip2 * total[pos];
      } else if (lcr == 'r') {
        ddos1[iado] += cr * total[pos] * d->dipole.pdip2;
      } else if (lcr == 'c') {
        ddos1[iado] += cl * d->dipole.pdip2 * total[pos] -
                       cr * total[pos] * d->dipole.pdip2;
      }
    }
  }

  {
    pos = d->g_index_list_pos[iado](5) - 1;
    if (pos != -1) {
      double n1 = key(0);
      double n2 = key(1);
      complex<double> sn =
          2. * sqrt(n1 * n2) / d->coef_abs[0] / d->coef_abs[1];
      complex<double> cl = d->dipole.bdip2[0] * d->alp2 * sn *
                           d->coef_lft[0] * d->coef_lft[1];
      complex<double> cr = d->dipole.bdip2[0] * d->alp2 * sn *
                           d->coef_rht[0] * d->coef_rht[1];
      if (lcr == 'l') {
        ddos1[iado] += cl * d->dipole.pdip2 * total[pos];
      } else if (lcr == 'r') {
        ddos1[iado] += cr * total[pos] * d->dipole.pdip2;
      } else if (lcr == 'c') {
        ddos1[iado] += cl * d->dipole.pdip2 * total[pos] -
                       cr * total[pos] * d->dipole.pdip2;
      }
    }
  }

  for (int mp = 0; mp < 2; ++mp) {
    pos = d->g_index_list_pos[iado](6 + mp) - 1;
    if (pos != -1) {
      double n1 = key(mp);
      complex<double> sn =
          d->dipole.bdip1[0] * d->alp1 * d->coef_abs[mp] * sqrt(n1 + 1);
      if (lcr == 'l') {
        ddos1[iado] += sn * (d->dipole.pdip1 * total[pos]);
      } else if (lcr == 'r') {
        ddos1[iado] += sn * (total[pos] * d->dipole.pdip1);
      } else if (lcr == 'c') {
        ddos1[iado] += sn * (d->dipole.pdip1 * total[pos] -
                             total[pos] * d->dipole.pdip1);
      }
    }
  }

  for (int mp = 0; mp < 2; ++mp) {
    pos = d->g_index_list_pos[iado](8 + mp) - 1;
    if (pos != -1) {
      double n1 = key(mp);
      complex<double> sn =
          d->dipole.bdip1[0] * d->alp1 * sqrt(n1) / d->coef_abs[mp];
      complex<double> cl = sn * d->coef_lft[mp];
      complex<double> cr = sn * d->coef_rht[mp];
      if (lcr == 'l') {
        ddos1[iado] += cl * d->dipole.pdip1 * total[pos];
      } else if (lcr == 'r') {
        ddos1[iado] += cr * total[pos] * d->dipole.pdip1;
      } else if (lcr == 'c') {
        ddos1[iado] += cl * d->dipole.pdip1 * total[pos] -
                       cr * total[pos] * d->dipole.pdip1;
      }
    }
  }

  for (int mp = 0; mp < 2; ++mp) {
    pos = d->g_index_list_pos[iado](10 + mp) - 1;
    if (pos != -1) {
      double n1 = key(mp);
      complex<double> sn = d->dipole.bdip2[0] * d->alp2 *
                           sqrt((n1 + 2) * (n1 + 1)) * d->coef_abs[mp] *
                           d->coef_abs[mp];
      if (lcr == 'l') {
        ddos1[iado] += sn * (d->dipole.pdip2 * total[pos]);
      } else if (lcr == 'r') {
        ddos1[iado] += sn * (total[pos] * d->dipole.pdip2);
      } else if (lcr == 'c') {
        ddos1[iado] += sn * (d->dipole.pdip2 * total[pos] -
                             total[pos] * d->dipole.pdip2);
      }
    }
  }

  {
    pos = d->g_index_list_pos[iado](12) - 1;
    if (pos != -1) {
      double n1 = key(0);
      double n2 = key(1);
      complex<double> sn = d->dipole.bdip2[0] * d->alp2 * 2.0 *
                           sqrt((n1 + 1) * (n2 + 1)) * d->coef_abs[0] *
                           d->coef_abs[1];
      if (lcr == 'l') {
        ddos1[iado] += sn * (d->dipole.pdip2 * total[pos]);
      } else if (lcr == 'r') {
        ddos1[iado] += sn * (total[pos] * d->dipole.pdip2);
      } else if (lcr == 'c') {
        ddos1[iado] += sn * (d->dipole.pdip2 * total[pos] -
                             total[pos] * d->dipole.pdip2);
      }
    }
  }
}

void rem_cal_filter_p(nNNmat& ddos1, const nNNmat& total,
                      const DEOM_DATA* d, const int iado) {
  key_vec key(d->keys[iado]);
  complex<double> gams = 0.0;
  for (int i = 2; i < d->nind; i++) gams += (double)(key(i)) * d->expn(i);
  int pos = d->g_index_list_pos[iado](0) - 1;

  ddos1[iado].setZero();

  if (pos != -1) {
    ddos1[iado] += -gams * total[pos];
    ddos1[iado] -= d->i * (d->hamt * total[pos] - total[pos] * d->hamt);

    double n1 = key(0);
    double n2 = key(1);
    complex<double> cl = -d->i * ((2 * n1 + 1) * d->coef_lft[0] +
                                  (2 * n2 + 1) * d->coef_lft[1]);
    complex<double> cr = -d->i * ((2 * n1 + 1) * d->coef_rht[0] +
                                  (2 * n2 + 1) * d->coef_rht[1]);
    ddos1[iado] += d->coef0 * (n1 - n2) * total[pos];
    ddos1[iado] += d->coef3 * (cl * total[pos] - cr * total[pos]);
    ddos1[iado] -=
        d->i * d->alp0 * (d->qmdt * total[pos] - total[pos] * d->qmdt);
    ddos1[iado] +=
        d->alp2 * (cl * d->qmdt * total[pos] - cr * total[pos] * d->qmdt);
  }

  pos = d->g_index_list_pos[iado](1) - 1;
  if (pos != -1) {
    double n1 = key(0);
    double n2 = key(1);
    complex<double> sn =
        sqrt(n1 * (n2 + 1)) / d->coef_abs[0] * d->coef_abs[1];
    complex<double> cl = -d->i * 2. * sn * d->coef_lft[0];
    complex<double> cr = -d->i * 2. * sn * d->coef_rht[0];
    ddos1[iado] += d->coef1 * sn * total[pos];
    ddos1[iado] += d->coef3 * (cl * total[pos] - cr * total[pos]);
    ddos1[iado] +=
        d->alp2 * (cl * d->qmdt * total[pos] - cr * total[pos] * d->qmdt);
  }

  pos = d->g_index_list_pos[iado](2) - 1;
  if (pos != -1) {
    double n1 = key(0);
    double n2 = key(1);
    complex<double> sn =
        sqrt((n1 + 1) * n2) * d->coef_abs[0] / d->coef_abs[1];
    complex<double> cl = -d->i * 2. * sn * d->coef_lft[1];
    complex<double> cr = -d->i * 2. * sn * d->coef_rht[1];
    ddos1[iado] += d->coef2 * sn * total[pos];
    ddos1[iado] += d->coef3 * (cl * total[pos] - cr * total[pos]);
    ddos1[iado] +=
        d->alp2 * (cl * d->qmdt * total[pos] - cr * total[pos] * d->qmdt);
  }

  for (int mp = 0; mp < 2; ++mp) {
    pos = d->g_index_list_pos[iado](mp + 3) - 1;
    if (pos != -1) {
      double n1 = key(mp);
      complex<double> sn =
          -d->i * sqrt(n1 * (n1 - 1)) / d->coef_abs[mp] / d->coef_abs[mp];
      complex<double> cl = sn * d->coef_lft[mp] * d->coef_lft[mp];
      complex<double> cr = sn * d->coef_rht[mp] * d->coef_rht[mp];
      ddos1[iado] += d->coef3 * (cl * total[pos] - cr * total[pos]);
      ddos1[iado] += d->alp2 * (cl * d->qmdt * total[pos] -
                                cr * total[pos] * d->qmdt);
    }
  }

  pos = d->g_index_list_pos[iado](5) - 1;
  if (pos != -1) {
    double n1 = key(0);
    double n2 = key(1);
    complex<double> sn = -d->i * 2. * sqrt(double(n1 * n2)) /
                         d->coef_abs[0] / d->coef_abs[1];
    const complex<double> cl = sn * d->coef_lft[0] * d->coef_lft[1];
    const complex<double> cr = sn * d->coef_rht[0] * d->coef_rht[1];
    ddos1[iado] += d->coef3 * (cl * total[pos] - cr * total[pos]);
    ddos1[iado] +=
        d->alp2 * (cl * d->qmdt * total[pos] - cr * total[pos] * d->qmdt);
  }

  for (int mp = 0; mp < 2; ++mp) {
    pos = d->g_index_list_pos[iado](6 + mp) - 1;
    if (pos != -1) {
      double n1 = key(mp);
      complex<double> sn =
          -d->i * d->alp1 * d->coef_abs[mp] * sqrt((double)(n1 + 1));
      ddos1[iado] += sn * (d->qmdt * total[pos] - total[pos] * d->qmdt);
    }
  }

  for (int mp = 0; mp < 2; ++mp) {
    pos = d->g_index_list_pos[iado](8 + mp) - 1;
    if (pos != -1) {
      double n1 = key(mp);
      complex<double> sn = -d->i * d->alp1 * sqrt((n1)) / d->coef_abs[mp];
      complex<double> cl = sn * d->coef_lft[mp];
      complex<double> cr = sn * d->coef_rht[mp];
      ddos1[iado] += cl * d->qmdt * total[pos] - cr * total[pos] * d->qmdt;
    }
  }

  for (int mp = 0; mp < 2; ++mp) {
    pos = d->g_index_list_pos[iado](10 + mp) - 1;
    if (pos != -1) {
      double n1 = key(mp);
      complex<double> sn = -d->i * d->alp2 * sqrt(((n1 + 2) * (n1 + 1))) *
                           d->coef_abs[mp] * d->coef_abs[mp];
      ddos1[iado] += sn * (d->qmdt * total[pos] - total[pos] * d->qmdt);
    }
  }

  pos = d->g_index_list_pos[iado](12) - 1;
  if (pos != -1) {
    double n1 = key(0);
    double n2 = key(1);
    complex<double> sn = -d->i * d->alp2 * 2.0 *
                         sqrt(((n1 + 1) * (n2 + 1))) * d->coef_abs[0] *
                         d->coef_abs[1];
    ddos1[iado] += sn * (d->qmdt * total[pos] - total[pos] * d->qmdt);
  }

  for (int mp1 = 0; mp1 < 2; ++mp1) {
    for (int mp2 = 2; mp2 < d->nind; ++mp2) {
      pos = d->g_index_list_pos[iado](11 + mp1 * (d->nind - 2) + mp2) - 1;
      if (pos != -1) {
        double n1 = key(mp1);
        double n2 = key(mp2);
        complex<double> sn =
            -d->i * (d->coef_lft[mp1] - d->coef_rht[mp1]) *
            sqrt((n1 * (n2 + 1))) / d->coef_abs[mp1] * d->coef_abs[mp2];
        ddos1[iado] += sn * total[pos];
      }
    }
  }

  for (int mp1 = 0; mp1 < 2; ++mp1) {
    for (int mp2 = 2; mp2 < d->nind; ++mp2) {
      pos =
          d->g_index_list_pos[iado](11 + (mp1 + 2) * (d->nind - 2) + mp2) -
          1;
      if (pos != -1) {
        double n1 = key[mp1];
        double n2 = key[mp2];
        complex<double> sn =
            -d->i * (d->coef_lft[mp2] - d->coef_rht[mp2]) *
            sqrt((n1 + 1) * n2) * d->coef_abs[mp1] / d->coef_abs[mp2];
        ddos1[iado] += sn * total[pos];
      }
    }
  }

  for (int mp1 = 0; mp1 < 2; ++mp1) {
    for (int mp2 = 2; mp2 < d->nind; ++mp2) {
      pos =
          d->g_index_list_pos[iado](11 + (mp1 + 4) * (d->nind - 2) + mp2) -
          1;
      if (pos != -1) {
        double n1 = key[mp1];
        double n2 = key[mp2];
        complex<double> sn = -d->i *
                             (d->coef_lft[mp1] * d->coef_lft[mp2] -
                              d->coef_rht[mp1] * d->coef_rht[mp2]) *
                             sqrt(n1 * n2) / d->coef_abs[mp2] /
                             d->coef_abs[mp1];
        ddos1[iado] += sn * total[pos];
      }
    }
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
                           const DEOM_SYL* syl, const int iado) {}

void syl_cal_filter_p(nNNmat& ddos1, const nNNmat& total,
                      const DEOM_DATA* d, const DEOM_SYL* syl,
                      const int iado) {}

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
  sum_Matrix.setZero();
  int pos = tree->find(1);

  if (pos != -1) {
    for (int nmod = 0; nmod < d->nmod; nmod++) {
      sum_Matrix += d->dipole1.sdip[nmod] * d->ddos[pos];
    }

    complex<double> cl =
        d->dipole1.bdip2[0] * d->alp2 * (d->coef_lft[0] + d->coef_lft[1]);
    sum_Matrix +=
        d->alp0 * d->dipole1.bdip0[0] * d->dipole1.pdip0 * d->ddos[pos] +
        cl * d->dipole1.pdip2 * d->ddos[pos];
  }

  for (int mp = 0; mp < 2; ++mp) {
    int hash_val = generate_key_plus(d->zerokey, d->emptykey, mp, d);
    pos = tree->find(hash_val);
    if (pos != -1) {
      complex<double> sn = d->dipole1.bdip1[0] * d->alp1 * d->coef_abs[mp];
      sum_Matrix += sn * d->dipole1.pdip1 * d->ddos[pos];
    }
  }

  for (int mp = 0; mp < 2; ++mp) {
    generate_key_plus(d->zerokey, d->emptykey, mp, d);
    int hash_val = generate_hash_value_plus(d->emptykey, mp, d);
    pos = tree->find(hash_val);
    if (pos != -1) {
      complex<double> sn = d->dipole1.bdip2[0] * d->alp2 *
                           sqrt((double)(2)) * d->coef_abs[mp] *
                           d->coef_abs[mp];
      sum_Matrix += sn * d->dipole1.pdip2 * d->ddos[pos];
    }
  }

  generate_key_plus(d->zerokey, d->emptykey, 0, d);
  int hash_val = generate_hash_value_plus(d->emptykey, 1, d);
  pos = tree->find(hash_val);
  if (pos != -1) {
    complex<double> sn = d->dipole1.bdip2[0] * d->alp2 * 2.0 *
                         d->coef_abs[0] * d->coef_abs[1];
    sum_Matrix += sn * d->dipole1.pdip2 * d->ddos[pos];
  }

  return sum_Matrix.trace();
}