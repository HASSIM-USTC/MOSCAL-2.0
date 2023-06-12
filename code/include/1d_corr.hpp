#pragma once
#include "deom.hpp"
#include "mr_level_1.hpp"
#include "print.hpp"

void od_corr(DEOM_DATA *d, const CTRL *c) {
  printf("task: cal 1d corr <\\mu(t)\\mu_1(0)>, only system operater is allowed!\n");

  int i = 0;
  oprt(c->lcr, i, d, rem_oprt);

  for (int ii = 0; ii < c->nt; ii++) {
    complex<double> trace = caculate_trace(d);
    print(c->fpolar, trace, c->ti + ii * c->dt);
    runge_kutta_4(i, d, c, rem_cal_1);
  }
}

void od_corr_filter(DEOM_DATA *d, const CTRL *c) {
  printf("task: cal 1d corr <\\mu(t)\\mu_1(0)>, dt method, only system operater is allowed!\n");

  int i = 0;
  oprt_filter(c->lcr, i, d, rem_oprt_filter);

  for (int ii = 0; ii < c->nt; ii++) {
    complex<double> trace = caculate_trace_filter(d);
    print(c->fpolar, trace, c->ti + ii * c->dt);
    init_nddo_clear_tree(d);
    runge_kutta_4_filter(i, d, c, rem_cal_1_filter);
  }
}

void od_corr_hei_filter(DEOM_DATA *d, const CTRL *c) {
  printf("task: cal 1d corr <\\mu(t)\\mu_1(0)>, only system operater is allowed!\n");

  int i = 0;
  allocator_system(i, d, c->lcr);
  Init_ddos_hei(d, d->ddos);

  for (int ii = 0; ii < c->nt; ii++) {
    complex<double> trace = caculate_trace_hei_filter(d);
    print(c->fpolar, trace, c->ti + ii * c->dt);
    init_nddo_clear_tree(d);
    runge_kutta_4_filter(i, d, c, rem_cal_1_hei_filter);
  }
}

void od_corr_sc2(DEOM_DATA *d, const CTRL *c, DEOM_SYL *syl) {
  printf("task: cal 1d corr <\\mu(t)\\mu_1(0)>.\n");

  int i = 0;
  gen_syl(d, syl);
  gen_keys(d);
  oprt(c->lcr, i, d, rem_oprt);
  A_assignments_B(i, d->nddo, d->ddos2, d->ddos1);

  for (int wii = 0; wii < c->w_len; wii++) {
    renew_check_flag(d);
    generate_syl_omega(d, syl, c->w_l[wii]);
    while (d->if_converged.flag) {
      check_if_converged(c, d);
      sc2_ddos_cal(d, syl, i, syl_cal_1);
    }
    print(c->fpolar, d->if_converged.trace, c->w_l[wii]);
  }
  generate_syl_omega(d, syl, 0);
}

void od_corr_sc2_filter(DEOM_DATA *d, const CTRL *c, DEOM_SYL *syl) {
  printf("task: cal 1d corr <\\mu(t)\\mu_1(0)>, filter.\n");

  int i = 0;
  gen_syl(d, syl);
  oprt_filter(c->lcr, i, d, rem_oprt_filter);

  for (int wii = 0; wii < c->w_len; wii++) {
    renew_check_flag(d);
    generate_syl_omega(d, syl, c->w_l[wii]);
    while (d->if_converged.flag) {
      init_nddo_clear_tree(d);
      sci_calculator_equal_ddos_filter_1(d, i, filter_2_p);
      check_if_converged_filter(c, d);
      sci_calculator_equal_ddos_filter_2(d, syl, i, syl_cal_1_filter);
    }
    print(c->fpolar, d->if_converged.trace, c->w_l[wii]);
  }
  generate_syl_omega(d, syl, 0);
}

void od_corr_Hei_sc2_filter(DEOM_DATA *d, const CTRL *c, DEOM_SYL *syl) {
  printf("task: cal 1d corr <\\mu(t)\\mu_1(0)>, filter, heisenberg picture.\n");

  int i = 0;
  allocator_system(i, d, c->lcr);
  gen_syl(d, syl);
  Init_ddos_hei(d, d->ddos2);

  for (int wii = 0; wii < c->w_len; wii++) {
    renew_check_flag(d);
    generate_syl_omega(d, syl, c->w_l[wii]);
    while (d->if_converged.flag) {
      init_nddo_clear_tree(d);
      sci_calculator_equal_ddos_filter_1(d, i, filter_hei_p);
      check_if_converged_hei_filter(c, d);
      sci_calculator_equal_ddos_filter_2(d, syl, i, syl_cal_hei_filter);
    }
    print(c->fpolar, d->if_converged.trace, c->w_l[wii]);
  }
  generate_syl_omega(d, syl, 0);
}

void od_corr_noise(DEOM_DATA *d, const CTRL *c) {
  printf("task: cal 1d corr <\\mu(t)\\mu_1(0)>, only bath operater is allowed!\n");

  int i = 0;
  oprt(c->lcr, i, d, rem_oprt);

  for (int ii = 0; ii < c->nt; ii++) {
    complex<double> trace = caculate_trace(d);
    print(c->fpolar, trace, c->ti + ii * c->dt);
    runge_kutta_4(i, d, c, rem_cal);
  }
}

void od_corr_noise_filter(DEOM_DATA *d, const CTRL *c) {
  printf("task: cal 1d corr <\\mu(t)\\mu_1(0)>, filter, only bath operater is allowed!\n");

  int i = 0;
  oprt_filter(c->lcr, i, d, rem_oprt_filter);

  for (int ii = 0; ii < c->nt; ii++) {
    complex<double> trace = caculate_trace_filter(d);
    print(c->fpolar, trace, c->ti + ii * c->dt);
    init_nddo_clear_tree(d);
    runge_kutta_4_filter(i, d, c, rem_cal_filter);
  }
}

void od_corr_noise_hei_filter(DEOM_DATA *d, const CTRL *c) {
  printf("task: cal 1d corr <\\mu(t)\\mu_1(0)>, only system operater is allowed!\n");

  int i = 0;
  allocator_bath(i, d, c->lcr);
  Init_ddos_hei(d, d->ddos);

  for (int ii = 0; ii < c->nt; ii++) {
    complex<double> trace = caculate_trace_hei_filter(d);
    print(c->fpolar, trace, c->ti + ii * c->dt);
    init_nddo_clear_tree(d);
    runge_kutta_4_filter(i, d, c, rem_cal_hei_filter);
  }
}

void od_corr_sc2_noise(DEOM_DATA *d, const CTRL *c, DEOM_SYL *syl) {
  printf("task: cal 1d corr <\\mu(t)\\mu_1(0)>, only bath operater is allowed!\n");

  int i = 0;
  gen_syl(d, syl);
  gen_keys(d);
  oprt(c->lcr, i, d, rem_oprt);
  A_assignments_B(i, d->nddo, d->ddos2, d->ddos1);

  for (int wii = 0; wii < c->w_len; wii++) {
    renew_check_flag(d);
    generate_syl_omega(d, syl, c->w_l[wii]);
    while (d->if_converged.flag) {
      check_if_converged(c, d);
      sc2_ddos_cal(d, syl, i, syl_ddos_cal);
    }
    print(c->fpolar, d->if_converged.trace, c->w_l[wii]);
  }
  generate_syl_omega(d, syl, 0);
}

void od_corr_sc2_noise_filter(DEOM_DATA *d, const CTRL *c, DEOM_SYL *syl) {
  printf("task: cal 1d corr <\\mu(t)\\mu_1(0)>, filter, only bath operater is allowed!\n");

  int i = 0;
  gen_syl(d, syl);
  oprt_filter(c->lcr, i, d, rem_oprt_filter);

  for (int wii = 0; wii < c->w_len; wii++) {
    renew_check_flag(d);
    generate_syl_omega(d, syl, c->w_l[wii]);
    while (d->if_converged.flag) {
      init_nddo_clear_tree(d);
      sci_calculator_equal_ddos_filter_1(d, i, filter_2_p);
      check_if_converged_filter(c, d);
      sci_calculator_equal_ddos_filter_2(d, syl, i, syl_cal_ddos_filter);
    }
    print(c->fpolar, d->if_converged.trace, c->w_l[wii]);
  }
  generate_syl_omega(d, syl, 0);
}

void od_corr_Hei_sc2_noise_filter(DEOM_DATA *d, const CTRL *c, DEOM_SYL *syl) {
  printf("task: cal 1d corr <\\mu(t)\\mu_1(0)>, filter, heisenberg picture, only system operater is allowed!\n");

  int i = 0;
  allocator_bath(i, d, c->lcr);
  Init_ddos_hei(d, d->ddos2);
  gen_syl(d, syl);

  for (int wii = 0; wii < c->w_len; wii++) {
    renew_check_flag(d);
    generate_syl_omega(d, syl, c->w_l[wii]);
    while (d->if_converged.flag) {
      init_nddo_clear_tree(d);
      sci_calculator_equal_ddos_filter_1(d, i, filter_hei_p);
      check_if_converged_hei_filter(c, d);
      sci_calculator_equal_ddos_filter_2(d, syl, i, syl_cal_1_hei_filter);
    }
    print(c->fpolar, d->if_converged.trace, c->w_l[wii]);
  }
  generate_syl_omega(d, syl, 0);
}
