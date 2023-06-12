#pragma once
#include "deom.hpp"

void Init_syl(DEOM_SYL* syl, DEOM_DATA* d, const Json& json) {
  syl->OMG = json["OMG"].get<double>();
  syl->nind = json["nind"].get<double>();
  syl->lwsg = json["lwsg"].get<double>();

  ComplexEigenSolver<MatrixXcd> eigensolver(
      (MatrixXcd)(d->i * d->hamt - syl->OMG * d->I));
  syl->S = eigensolver.eigenvalues();
  MatrixXcd V = eigensolver.eigenvectors();
  MatrixXcd VI = V.inverse();
  allocate_syl(syl->V, syl->VI, V, VI);

  syl->twsg.resize(syl->nind, syl->lwsg);
  for (int i_nind = 0; i_nind < syl->nind; i_nind++)
    for (int i = 0; i < syl->lwsg; i++)
      syl->twsg(i_nind, i) =
          json["twsg"][INDEX2(i_nind, i, syl->lwsg)].get<double>();

  syl->expn.resize(syl->nind);
  for (int i_nind = 0; i_nind < syl->nind; i_nind++) {
    syl->expn(i_nind) = d->expn(syl->twsg(i_nind, 0));
    for (int i = 1; i < syl->lwsg; i++) {
      if (abs(syl->expn(i_nind) - d->expn(syl->twsg(i_nind, i))) > 1e-7) {
        printf("invaild syl input, exit\n");
        exit(41);
      }
    }
  }
}

void Init_syl_omg(DEOM_SYL* syl, DEOM_DATA* d, const Json& json) {
  syl->OMG = json["OMG"].get<double>();
  int syl_size = (int)d->comb_list(d->lmax + syl->nind, syl->nind);

  syl->R.resize(syl_size);
  syl->U.resize(syl_size);
  syl->UI.resize(syl_size);
  syl->keys.resize(syl_size, syl->nind);
  syl->keys.setZero();

  for (int iado_ = 0; iado_ < syl_size; iado_++) {
    syl->R[iado_].resize(NSYS);
    syl->U[iado_].resize(NSYS, NSYS);
    syl->UI[iado_].resize(NSYS, NSYS);
  }
}

complex<double> element(const Json& json) {
  return json["real"].get<double>() +
         complex<double>(0, 1) * json["imag"].get<double>();
}

complex<double> element(int index, const Json& json) {
  return json["real"][index].get<double>() +
         complex<double>(0, 1) * json["imag"][index].get<double>();
}

void json2vector(VectorXcd& vec, int size_vector, const Json& json) {
  vec = VectorXcd::Zero(size_vector);
  if (json["if_initial"].get<bool>()) {
    for (int i = 0; i < size_vector; i++) {
      vec(i) = element(i, json);
    }
  }
}

void json2matrix(MatrixNcd& mat, const Json& json) {
  mat.resize(NSYS, NSYS);
  if (json["if_initial"].get<bool>()) {
    for (int i = 0; i < NSYS; i++) {
      for (int j = 0; j < NSYS; j++) {
        insert(mat, i, j, element(INDEX2(i, j), json));
      }
    }
  } else {
    mat.setZero();
  }
  init_nNNmat(mat);
}

void json2matrices(nNNmat& matrix, int size_array, const Json& json) {
  matrix = nNNmat(size_array);
  for (int k = 0; k < size_array; k++) {
    matrix[k].resize(NSYS, NSYS);
    if (json["if_initial"].get<bool>()) {
      for (int i = 0; i < NSYS; i++)
        for (int j = 0; j < NSYS; j++)
          insert(matrix[k], i, j, element(INDEX3(k, i, j), json));
      init_nNNmat(matrix[k]);
    } else {
      matrix[k].setZero();
    }
  }
}

template <typename T>
void allocate_qmdt(T* allocate, DEOM_DATA* d, const Json& json) {
  json2matrices(allocate->qmdta_l, d->nind,
                json.value("qmdta_l", d->empty_list));
  json2matrices(allocate->qmdta_r, d->nind,
                json.value("qmdta_r", d->empty_list));
  json2matrices(allocate->qmdtc_l, d->nind,
                json.value("qmdtc_l", d->empty_list));
  json2matrices(allocate->qmdtc_r, d->nind,
                json.value("qmdtc_r", d->empty_list));
  json2matrices(allocate->qmdt2a_l, d->nind * d->nind,
                json.value("qmdt2a_l", d->empty_list));
  json2matrices(allocate->qmdt2a_r, d->nind * d->nind,
                json.value("qmdt2a_r", d->empty_list));
  json2matrices(allocate->qmdt2b_l, d->nind * d->nind,
                json.value("qmdt2b_l", d->empty_list));
  json2matrices(allocate->qmdt2b_r, d->nind * d->nind,
                json.value("qmdt2b_r", d->empty_list));
  json2matrices(allocate->qmdt2c_l, d->nind * d->nind,
                json.value("qmdt2c_l", d->empty_list));
  json2matrices(allocate->qmdt2c_r, d->nind * d->nind,
                json.value("qmdt2c_r", d->empty_list));
  json2matrices(allocate->renormalize, d->nmod,
                json.value("renormalize", d->empty_list));
}

void Init_dip(DIPOLE* dipole, DEOM_DATA* d, const Json& json) {
  dipole->on = true;
  json2matrices(dipole->sdip, d->nmod,
                json.value("sdip_cub", d->empty_list));
  json2vector(dipole->bdip0, d->nind,
              json.value("bdip0_cub", d->empty_list));
  json2vector(dipole->bdip1, d->nind,
              json.value("bdip1_cub", d->empty_list));
  json2vector(dipole->bdip2, d->nind,
              json.value("bdip2_cub", d->empty_list));
  json2matrix(dipole->pdip0, json.value("pdip0_cub", d->empty_list));
  json2matrix(dipole->pdip1, json.value("pdip1_cub", d->empty_list));
  json2matrix(dipole->pdip2, json.value("pdip2_cub", d->empty_list));
  allocate_qmdt(dipole, d, json);
}

void Init_sys(DEOM_DATA* d, const Json& json) {
  printf("MOSCAL 2.1\n");
  printf("Author: Zi-Hao Chen : czh5@mail.ustc.edu.cn\n");
  printf("        Yu Su       : suyupilemao@mail.ustc.edu.cn\n");
  printf("        Yao Wang    : wy2010@ustc.edu.cn\n");
  printf("        YiJing Yan  : yjyan@ustc.edu.cn\n");
  printf(
      "Update: some of the init part are move to python script. So the "
      "input etar/etal must be normalization and the print function has "
      "been banned.\n");

#if defined(STD)
  printf("This part of code not support parallel filter.\n");
#endif  // STD

  d->empty_list["if_initial"] = false;
  d->empty_element["real"] = 0;
  d->empty_element["imag"] = 0;

  d->nmax = json["nmax"].get<int>();
  d->nind = json["nind"].get<int>();
  d->lmax = json["lmax"].get<int>();
  d->lmax_fp = json.value("lmax_fp", 0);
  d->lmax_ma = json.value("lmax_ma", 0);
  d->nmod = json["nmod"].get<int>();
  d->ferr = json.value("ferr", 1e-10);
  d->alp0 = json.value("alp0", 0.0);
  d->alp1 = json.value("alp1", 0.0);
  d->alp2 = json.value("alp2", 0.0);
  d->length_alp = json.value("length_alp", 0.0);
  json2vector(d->alp_list, d->length_alp,
              json.value("alp_list", d->empty_list));
  //
  d->coef0 = element(json.value("coef0", d->empty_element));
  d->coef1 = element(json.value("coef1", d->empty_element));
  d->coef2 = element(json.value("coef2", d->empty_element));
  d->coef3 = element(json.value("coef3", d->empty_element));
  d->omgb_BO = element(json.value("omgb_BO", d->empty_element));
  d->one_over_2_times_mass =
      element(json.value("one_over_2_times_mass", d->empty_element));

  int combmax = d->nind + d->lmax + 4;
  d->comb_list.resize(combmax, combmax);
  for (int j = 0; j < combmax; j++) {
    d->comb_list(0, j) = 0;
  }
  d->comb_list(0, 0) = 1;

  for (int i = 1; i < combmax; i++) {
    for (int j = 1; j < combmax; j++) {
      d->comb_list(i, j) =
          d->comb_list((i - 1), j) + d->comb_list((i - 1), j - 1);
    }
    d->comb_list(i, 0) = 1;
  }

  if (d->comb_list((d->lmax + d->nind), d->lmax) < d->nmax) {
    printf("I have set nmax to %lld, not %d.\n",
           (d->comb_list(d->lmax + d->nind, d->lmax)) +
               (d->comb_list(d->lmax + d->nind, d->lmax)) / 5,
           d->nmax);
    d->nmax = d->comb_list(d->lmax + d->nind, d->lmax) +
              (d->comb_list(d->lmax + d->nind, d->lmax)) / 5;
  }

  d->emptykey = key_vec(d->nind);
  d->zerokey = key_vec(d->nind);
  d->tempkey = key_vec(d->nind);
  d->emptykey.setZero();
  d->zerokey.setZero();
  d->tempkey.setZero();

  d->keys = array_key_vec(d->nmax);
  for (int iado = 0; iado < d->nmax; iado++) {
    d->keys[iado].resize(d->nind);
    d->keys[iado].setZero();
  }

  d->ddos = nNNmat(d->nmax);
  for (int iado_ = 0; iado_ < d->nmax; iado_++) {
    d->ddos[iado_].resize(NSYS, NSYS);
    d->ddos[iado_].setZero();
  }

  d->tree_hei = new Trie(d->nmax);
  d->tree = new Trie(d->nmax);

  if (json.value("read_rho0", false)) {
    json2matrix(d->ddos[0], json["rho0"]);
  } else {
    int inistate = json.value("inistate", 0);
    insert(d->ddos[0], inistate, inistate, 1.0);
  }

  d->nddo = 1;
  d->tree->try_insert(1, 0);
  printf("size of int: %lu\t size of llint :%lu\n", sizeof(int),
         sizeof(llint));

  json2matrix(d->ham1, json["ham1"]);
  json2matrix(d->qmd1, json.value("qmd1", d->empty_list));

  d->I.resize(NSYS, NSYS);
  d->I.setZero();
  for (int i = 0; i < NSYS; i++) {
    insert(d->I, i, i, 1);
  }

  json2vector(d->expn, d->nind, json["expn"]);
  json2vector(d->coef_lft, d->nind, json.value("coef_lft", d->empty_list));
  json2vector(d->coef_rht, d->nind, json.value("coef_rht", d->empty_list));
  json2vector(d->coef_abs, d->nind, json["coef_abs"]);

  for (int i = 0; i < d->nind; i++) {
    d->coef_abs[i] = sqrt(d->coef_abs[i]);
  }

  printf("qmdt are given by input\n");
  allocate_qmdt(d, d, json);

  d->qmdt2a_equal_0.resize(d->nind, d->nind);
  d->qmdt2b_equal_0.resize(d->nind, d->nind);
  d->qmdt2c_equal_0.resize(d->nind, d->nind);
  for (int mp1 = 0; mp1 < d->nind; mp1++) {
    for (int mp2 = 0; mp2 < d->nind; mp2++) {
      d->qmdt2a_equal_0(mp1, mp2) =
          is_valid(d->qmdt2a_l[INDEX2(mp1, mp2, d->nind)], 1e-15);
      d->qmdt2b_equal_0(mp1, mp2) =
          is_valid(d->qmdt2b_l[INDEX2(mp1, mp2, d->nind)], 1e-15);
      d->qmdt2c_equal_0(mp1, mp2) =
          is_valid(d->qmdt2c_l[INDEX2(mp1, mp2, d->nind)], 1e-15);
    }
  }

  d->hamt.resize(NSYS, NSYS);
  d->qmdt.resize(NSYS, NSYS);
  d->hamt = d->ham1;
  d->qmdt = d->qmd1;

  d->ham1_abs = d->ham1.squaredNorm();
  d->qmd1_abs.resize(d->nind);
  for (int i_nind = 0; i_nind < d->nind; i_nind++)
    d->qmd1_abs(i_nind) = d->qmdta_r[i_nind].squaredNorm();
}

void Init_cfw(DEOM_DATA* d, const Json& json) {
  d->cfw.tau = element(json["tau"]);
  d->cfw.alpha = element(json["alpha"]);
  d->cfw.tf = element(json["tf"]);
}

void Init_sys_filter(DEOM_DATA* d, const Json& json) {
  d->g_index_list_pos.resize(d->nmax);
  d->keys1 = array_key_vec(d->nmax);
  for (int iado = 0; iado < d->nmax; iado++) {
    d->keys1[iado].resize(d->nind);

#if defined(FERMI_QUAD)
#elif defined(FERMI_LINEAR)
    d->g_index_list_pos[iado].resize(d->nind + 1);
#elif defined(BOSE_QUAD)
    d->g_index_list_pos[iado].resize(2 * d->nind + 3 * d->nind * d->nind +
                                     1);
#elif defined(BOSE_LINEAR)
    d->g_index_list_pos[iado].resize(2 * d->nind + 1);
#elif defined(BSM)
    d->g_index_list_pos[iado].resize(6 * d->nind + 1);
#elif defined(BSM_ACTION)
    d->g_index_list_pos[iado].resize(6 * d->nind + 1);
#elif defined(MD)
    d->g_index_list_pos[iado].resize(6 * d->nind + 1);
#endif

    d->g_index_list_pos[iado].setZero();
    d->keys1[iado].setZero();
  }
}

void Init_bak(DEOM_BAK* b, DEOM_DATA* d) {
  b->ddos = nNNmat(d->nmax);
  for (int iado_ = 0; iado_ < d->nmax; iado_++) {
    b->ddos[iado_].resize(NSYS, NSYS);
    b->ddos[iado_] = d->ddos[iado_];
  }
  b->nddo = d->nddo;
}

void Init_bak(DEOM_DATA* d, DEOM_BAK* b) {
  for (int iado_ = 0; iado_ < d->nmax; iado_++) {
    d->ddos[iado_].setZero();
    d->ddos[iado_] = b->ddos[iado_];
  }
  d->nddo = b->nddo;
}

void Init_ctrl(CTRL* c, const Json& json) {
  c->dt = json.value("dt", 0.0);
  c->ti = json.value("ti", 0.0);
  c->tf = json.value("tf", 0.0);
  c->ferr = json.value("ferr", 0.0);
  c->ferr_adapt = json.value("ferr_adapt", 0.0);
  c->step = json.value("step", 0.0);
  c->lcr = json.value("lcr", "E").c_str()[0];

  if (c->ti < 0) {
    c->nt_i = ceil(abs(c->ti) / c->dt);
    c->nt_f = ceil(abs(c->tf) / c->dt);
    c->nt = c->nt_i + c->nt_f;
    c->ti = -c->nt_i * c->dt;
  } else {
    c->nt_i = ceil(abs(c->ti) / c->dt);
    c->nt_f = ceil(abs(c->tf) / c->dt);
    c->nt = c->nt_f - c->nt_i;
    c->ti = c->nt_i * c->dt;
  }

  c->w_len = json.value("w_len", 0);
  c->w_l.resize(c->w_len);
  for (int i = 0; i < c->w_len; i++)
    c->w_l(i) = json["w_l"][i].get<double>();
}

void Init_aux_spe_sc2(DEOM_DATA* d, const Json& json) {
  init_ddos(d->nmax, d->ddos1, d->ddos2);
}

void Init_aux_sc2(DEOM_DATA* d, const Json& json) {
  init_ddos(d->nmax, d->ddos1);
}

void Init_aux_dt(DEOM_DATA* d, const Json& json) {
  init_ddos(d->nmax, d->ddos1, d->ddos2, d->ddos3);
}

void Init_aux_dt_32(DEOM_DATA* d, const Json& json) {
  init_ddos(d->nmax, d->ddos1, d->ddos2, d->ddos3, d->ddos4);
}

void Init_aux_dt_54(DEOM_DATA* d, const Json& json) {
  init_ddos(d->nmax, d->ddos1, d->ddos2, d->ddos3, d->ddos4, d->ddos5,
            d->ddos6, d->ddos7, d->ddos8, d->ddos9, d->ddos10);
}
