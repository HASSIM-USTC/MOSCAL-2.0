#pragma once
#include "deom.hpp"
#include "mr_level_2.hpp"
#include "rk_include.hpp"

template <typename FUN>
void runge_kutta_4_filter(int i, DEOM_DATA *d, const CTRL *c,
                          FUN rem_cal_filter) {
#pragma omp parallel default(shared)
  {
    clean(i, d->nsave[6], d->g_index_list_pos, d->ddos1);
    filter(d->nsave[6], i, d->tree, d, filter_p);
    A_assignments_B(i, d->nsave[6], d->ddos, d->ddos1, d->keys, d->keys1);
    nsave_print(d->nsave, d->nddo, 1);

    construct_Mapping_filter(d->nsave[1], i, d->tree, d);
    nsave_print(d->nsave, d->nddo, 2);
    rem_cal_filter(d->ddos1, d->ddos, d->nsave[2], i, d);
    RK1_filter(i, d->nsave[1], d->nsave[2], d->ddos, d->ddos1, d->ddos3,
               c->dt / 2.0);

    construct_Mapping_filter(d->ddos3, d->nsave[2], i, d->tree, d);
    nsave_print(d->nsave, d->nddo, 3);
    rem_cal_filter(d->ddos2, d->ddos3, d->nsave[3], i, d);
    RK2_filter(i, d->nsave[1], d->nsave[2], d->nsave[3], d->ddos, d->ddos1,
               d->ddos2, d->ddos3, c->dt / 2.0);

    construct_Mapping_filter(d->ddos3, d->nsave[3], i, d->tree, d);
    nsave_print(d->nsave, d->nddo, 4);
    rem_cal_filter(d->ddos2, d->ddos3, d->nsave[4], i, d);
    RK3_filter(i, d->nsave[1], d->nsave[3], d->nsave[4], d->ddos, d->ddos1,
               d->ddos2, d->ddos3, c->dt);

    construct_Mapping_filter(d->ddos3, d->nsave[4], i, d->tree, d);
    nsave_print(d->nsave, d->nddo, 5);
    rem_cal_filter(d->ddos2, d->ddos3, d->nsave[5], i, d);
    RK4_filter(i, d->nsave[1], d->nsave[4], d->nsave[5], d->ddos, d->ddos1,
               d->ddos2, c->dt / 6.0);
  }
}

template <typename FUN>
void runge_kutta_4(int i, DEOM_DATA *d, const CTRL *c, FUN rem_cal) {
#pragma omp parallel default(shared)
  {
    rem_cal(d->ddos1, d->ddos, d->nddo, i, d);
    RK1(i, d->nddo, d->ddos, d->ddos1, d->ddos3, c->dt / 2.0);

    rem_cal(d->ddos2, d->ddos3, d->nddo, i, d);
    RK2(i, d->nddo, d->ddos, d->ddos1, d->ddos2, d->ddos3, c->dt / 2.0);

    rem_cal(d->ddos2, d->ddos3, d->nddo, i, d);
    RK3(i, d->nddo, d->ddos, d->ddos1, d->ddos2, d->ddos3, c->dt);

    rem_cal(d->ddos2, d->ddos3, d->nddo, i, d);
    RK4(i, d->nddo, d->ddos, d->ddos1, d->ddos2, d->ddos3, c->dt / 6.0);
  }
}

template <typename FUN1, typename FUN2>
void runge_kutta_4_filter_pulse(int i, int ii, DEOM_DATA *d, const CTRL *c,
                                FUN1 rem_cal_filter,
                                FUN2 pulse_change_system) {
#pragma omp parallel default(shared)
  {
    clean(i, d->nsave[6], d->g_index_list_pos, d->ddos1);
    filter(d->nsave[6], i, d->tree, d, filter_p);
    A_assignments_B(i, d->nsave[6], d->ddos, d->ddos1, d->keys, d->keys1);
    nsave_print(d->nsave, d->nddo, 1);

    construct_Mapping_filter(d->nsave[1], i, d->tree, d);
    nsave_print(d->nsave, d->nddo, 2);
    pulse_change_system(d, c->ti + ii * c->dt);
    rem_cal_filter(d->ddos1, d->ddos, d->nsave[2], i, d);
    RK1_filter(i, d->nsave[1], d->nsave[2], d->ddos, d->ddos1, d->ddos3,
               c->dt / 2.0);

    construct_Mapping_filter(d->ddos3, d->nsave[2], i, d->tree, d);
    nsave_print(d->nsave, d->nddo, 3);
    pulse_change_system(d, c->ti + ii * c->dt + c->dt / 2.0);
    rem_cal_filter(d->ddos2, d->ddos3, d->nsave[3], i, d);
    RK2_filter(i, d->nsave[1], d->nsave[2], d->nsave[3], d->ddos, d->ddos1,
               d->ddos2, d->ddos3, c->dt / 2.0);

    construct_Mapping_filter(d->ddos3, d->nsave[3], i, d->tree, d);
    nsave_print(d->nsave, d->nddo, 4);
    pulse_change_system(d, c->ti + ii * c->dt + c->dt / 2.0);
    rem_cal_filter(d->ddos2, d->ddos3, d->nsave[4], i, d);
    RK3_filter(i, d->nsave[1], d->nsave[3], d->nsave[4], d->ddos, d->ddos1,
               d->ddos2, d->ddos3, c->dt);

    construct_Mapping_filter(d->ddos3, d->nsave[4], i, d->tree, d);
    nsave_print(d->nsave, d->nddo, 5);
    pulse_change_system(d, c->ti + ii * c->dt + c->dt);
    rem_cal_filter(d->ddos2, d->ddos3, d->nsave[5], i, d);
    RK4_filter(i, d->nsave[1], d->nsave[4], d->nsave[5], d->ddos, d->ddos1,
               d->ddos2, c->dt / 6.0);
  }
}

template <typename FUN>
void runge_kutta_32_filter(int i, DEOM_DATA *d, CTRL *c,
                           FUN rem_cal_filter) {
#pragma omp parallel default(shared)
  {
    clean(i, d->nsave[6], d->g_index_list_pos, d->ddos1);
    filter(d->nsave[6], i, d->tree, d, filter_p);
    A_assignments_B(i, d->nsave[6], d->ddos, d->ddos1, d->keys, d->keys1);
    nsave_print(d->nsave, d->nddo, 1);

    construct_Mapping_filter(d->nsave[1], i, d->tree, d);
    nsave_print(d->nsave, d->nddo, 2);
    rem_cal_filter(d->ddos1, d->ddos, d->nsave[2], i, d);
    RK32_1_filter(i, d->nsave[1], d->nsave[2], d->ddos, d->ddos1, d->ddos3,
                  c->dt);

    construct_Mapping_filter(d->ddos3, d->nsave[2], i, d->tree, d);
    nsave_print(d->nsave, d->nddo, 3);
    rem_cal_filter(d->ddos2, d->ddos3, d->nsave[3], i, d);
    RK32_2_filter(i, d->nsave[1], d->nsave[2], d->nsave[3], d->ddos,
                  d->ddos1, d->ddos2, d->ddos3, d->ddos4, c->dt);

    construct_Mapping_filter(d->ddos3, d->nsave[3], i, d->tree, d);
    nsave_print(d->nsave, d->nddo, 4);
    rem_cal_filter(d->ddos2, d->ddos3, d->nsave[4], i, d);
    RK32_3_filter(i, d->nsave[1], d->nsave[3], d->nsave[4], d->ddos,
                  d->ddos1, d->ddos2, d->ddos3, d->ddos4, c->dt);

    construct_Mapping_filter(d->ddos3, d->nsave[4], i, d->tree, d);
    nsave_print(d->nsave, d->nddo, 5);
    rem_cal_filter(d->ddos2, d->ddos3, d->nsave[5], i, d);
    RK32_4_filter(i, d->nsave[4], d->nsave[5], d->ddos2, d->ddos3,
                  d->ddos4, c->dt);
  }
  adapt_rk(i, d->nsave[4], d->nsave[5], d->ddos, d->ddos3, d->ddos4, 0.33,
           c);
}

template <typename FUN>
void runge_kutta_54_filter(int i, DEOM_DATA *d, CTRL *c,
                           FUN rem_cal_filter) {
#pragma omp parallel default(shared)
  {
    clean(i, d->nsave[6], d->g_index_list_pos, d->ddos1);
    filter(d->nsave[6], i, d->tree, d, filter_p);
    A_assignments_B(i, d->nsave[6], d->ddos, d->ddos1, d->keys, d->keys1);
    nsave_print(d->nsave, d->nddo, 1);

    construct_Mapping_filter(d->nsave[1], i, d->tree, d);
    nsave_print(d->nsave, d->nddo, 2);
    rem_cal_filter(d->ddos1, d->ddos, d->nsave[2], i, d);
    RK54_1_filter(i, d->nsave, d->ddos, d->ddos1, d->ddos8, c->dt);

    construct_Mapping_filter(d->ddos8, d->nsave[2], i, d->tree, d);
    nsave_print(d->nsave, d->nddo, 3);
    rem_cal_filter(d->ddos2, d->ddos8, d->nsave[3], i, d);
    RK54_2_filter(i, d->nsave, d->ddos, d->ddos1, d->ddos2, d->ddos8,
                  d->ddos9, d->ddos10, c->dt);

    construct_Mapping_filter(d->ddos8, d->nsave[3], i, d->tree, d);
    nsave_print(d->nsave, d->nddo, 4);
    rem_cal_filter(d->ddos3, d->ddos8, d->nsave[4], i, d);
    RK54_3_filter(i, d->nsave, d->ddos, d->ddos1, d->ddos2, d->ddos3,
                  d->ddos8, d->ddos9, d->ddos10, c->dt);

    construct_Mapping_filter(d->ddos8, d->nsave[4], i, d->tree, d);
    nsave_print(d->nsave, d->nddo, 5);
    rem_cal_filter(d->ddos4, d->ddos8, d->nsave[5], i, d);
    RK54_4_filter(i, d->nsave, d->ddos, d->ddos1, d->ddos2, d->ddos3,
                  d->ddos4, d->ddos8, d->ddos9, d->ddos10, c->dt);

    construct_Mapping_filter(d->ddos8, d->nsave[5], i, d->tree, d);
    nsave_print(d->nsave, d->nddo, 6);
    rem_cal_filter(d->ddos5, d->ddos8, d->nsave[6], i, d);
    RK54_5_filter(i, d->nsave, d->ddos, d->ddos1, d->ddos2, d->ddos3,
                  d->ddos4, d->ddos5, d->ddos8, d->ddos9, d->ddos10,
                  c->dt);

    construct_Mapping_filter(d->ddos8, d->nsave[6], i, d->tree, d);
    nsave_print(d->nsave, d->nddo, 7);
    rem_cal_filter(d->ddos6, d->ddos8, d->nsave[7], i, d);
    RK54_6_filter(i, d->nsave, d->ddos6, d->ddos9, d->ddos10, c->dt);

    construct_Mapping_filter(d->ddos9, d->nsave[7], i, d->tree, d);
    nsave_print(d->nsave, d->nddo, 8);
    rem_cal_filter(d->ddos7, d->ddos9, d->nsave[8], i, d);
    RK54_7_filter(i, d->nsave, d->ddos7, d->ddos9, d->ddos10, c->dt);
  }
  adapt_rk(i, d->nsave[7], d->nsave[8], d->ddos, d->ddos9, d->ddos10, 0.2,
           c);
}