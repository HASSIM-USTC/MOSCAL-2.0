#pragma once
#include "deom.hpp"
#include "mr_level_1.hpp"

void equilibrium_md(DEOM_DATA *d, const CTRL *c) {
  printf("task: equilibrium, cal \\rho^{eq}, dt method, wtih filter\n");

  int i = 0;
  gen_keys(d);
  print_md(c->frho, c->fcur, 0, d, c);

  for (int ii = 0; ii < c->nt; ii++) {
    runge_kutta_4(i, d, c, rem_cal);
    print_md(c->frho, c->fcur, ii, d, c);
  }
}

void equilibrium(DEOM_DATA *d, const CTRL *c) {
  printf("task: equilibrium, cal \\rho^{eq}, dt-method\n");

  int i = 0;
  gen_keys(d);

  for (int ii = 0; ii < c->nt; ii++) {
    print(c->frho, c->fcur, ii, d, c);
    runge_kutta_4(i, d, c, rem_cal);
  }
}

void equilibrium_sc2(DEOM_DATA *d, const CTRL *c, DEOM_SYL *syl) {
  printf("task: equilibrium, cal \\rho^{eq}, sc2 method\n");

  int i = 0;
  gen_syl(d, syl);
  gen_keys(d);

  for (int ii = 0; ii < c->nt; ii++) {
    print(c->frho, c->fcur, ii, d, c);
    sc2_cal(d, syl, i, syl_cal);
  }
}

void equilibrium_md_filter(DEOM_DATA *d, const CTRL *c) {
  printf("task: equilibrium, cal \\rho^{eq}, dt method, wtih filter\n");

  int i = 0;
  print_md_filter(c->frho, c->fcur, 0, d, c);

  for (int ii = 0; ii < c->nt; ii++) {
    init_nddo_clear_tree(d);
    runge_kutta_4_filter(i, d, c, rem_cal_filter);
    print_md_filter(c->frho, c->fcur, ii, d, c);
  }
}

void equilibrium_filter(DEOM_DATA *d, const CTRL *c) {
  printf("task: equilibrium, cal \\rho^{eq}, dt method, wtih filter\n");

  int i = 0;
  print_filter(c->frho, c->fcur, 0, d, c);

  for (int ii = 0; ii < c->nt; ii++) {
    init_nddo_clear_tree(d);
    runge_kutta_4_filter(i, d, c, rem_cal_filter);
    print_filter(c->frho, c->fcur, ii, d, c);
  }
}

void equilibrium_32_filter(DEOM_DATA *d, CTRL *c) {
  printf("task: equilibrium, cal \\rho^{eq}, dt method, wtih filter\n");

  int i = 0;
  print_adapt_rk(c->frho, d, c);

  for (int ii = 0; ii < c->nt; ii++) {
    init_nddo_clear_tree(d);
    runge_kutta_32_filter(i, d, c, rem_cal_filter);
    print_adapt_rk(c->frho, d, c);
  }
}

void equilibrium_54_filter(DEOM_DATA *d, CTRL *c) {
  printf("task: equilibrium, cal \\rho^{eq}, dt method, wtih filter\n");

  int i = 0;
  print_adapt_rk(c->frho, d, c);

  for (int ii = 0; ii < c->nt; ii++) {
    init_nddo_clear_tree(d);
    runge_kutta_54_filter(i, d, c, rem_cal_filter);
    print_adapt_rk(c->frho, d, c);
  }
}

void equilibrium_sc2_filter(DEOM_DATA *d, const CTRL *c, DEOM_SYL *syl) {
  printf("task: equilibrium, cal \\rho^{eq}, sc2 method, wtih filter\n");

  int i = 0;
  gen_syl(d, syl);

  for (int ii = 0; ii < c->nt; ii++) {
    init_nddo_clear_tree(d);
    sc2_cal_filter_1(d, i, filter_p);
    sc2_cal_filter_2(d, syl, i, syl_cal_filter);
    print_filter(c->frho, c->fcur, ii, d, c);
  }
}
