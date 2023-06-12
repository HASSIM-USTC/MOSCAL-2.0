#pragma once
#include "deom.hpp"
#include "mr_level_1.hpp"

void equilibrium_imag_filter(DEOM_DATA *d, CTRL *c) {
  printf("task: imaginary time evolution, cal \\rho^{eq}\n");
  int i = 0;
  print_filter(c->frho, c->fcur, 0, d, c);

  for (int ii = 0; ii < c->nt; ii++) {
    init_nddo_clear_tree(d);
    runge_kutta_4_filter(i, d, c, rem_cal_imag_filter);
    print_filter(c->frho, c->fcur, ii, d, c);
  }
}

void equilibrium_cfw_forward_filter(DEOM_DATA *d, CTRL *c) {
  printf("task: cfw, forward, cal \\rho^{eq}\n");
  int i = 0;
  print_filter(c->frho, c->fcur, 0, d, c);

  for (int ii = 0; ii < c->nt; ii++) {
    init_nddo_clear_tree(d);
    runge_kutta_4_filter_pulse(i, ii, d, c, rem_cal_cfw_filter, lambda_forward);
    print_filter(c->frho, c->fcur, ii, d, c);
  }
}

void equilibrium_cfw_backward_filter(DEOM_DATA *d, CTRL *c) {
  printf("task: cfw, backward, cal \\rho^{eq}\n");
  int i = 0;
  print_filter(c->frho, c->fcur, 0, d, c);

  for (int ii = 0; ii < c->nt; ii++) {
    init_nddo_clear_tree(d);
    runge_kutta_4_filter_pulse(i, ii, d, c, rem_cal_cfw_filter,
                               lambda_backward);
    print_filter(c->frho, c->fcur, ii, d, c);
  }
}