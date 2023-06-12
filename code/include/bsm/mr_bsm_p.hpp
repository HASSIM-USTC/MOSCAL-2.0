#pragma once

#include "algebra.hpp"
#include "deom.hpp"
#include "linearalgebra.hpp"

void construct_Mapping_p(DEOM_DATA *d, int iado) {}

void rem_oprt_p(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d, const int iado,
              const char lcr) {}

void rem_cal_p(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
               const int iado) {}

void rem_cal_1_p(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
                 const int iado) {}

void syl_cal_p(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
               const DEOM_SYL *syl, const int iado) {}

void syl_cal_1_p(nNNmat &ddos1, const nNNmat &total, const nNNmat &ddos2,
                 const DEOM_DATA *d, const DEOM_SYL *syl, const int iado) {}

void syl_cal_ddos_p(nNNmat &ddos1, const nNNmat &total, const nNNmat &ddos2,
                    const DEOM_DATA *d, const DEOM_SYL *syl, const int iado) {}

complex<double> caculate_trace(DEOM_DATA *d) {
  MatrixXcd sum_Matrix = MatrixXcd::Zero(NSYS, NSYS);
  return sum_Matrix.trace();
}