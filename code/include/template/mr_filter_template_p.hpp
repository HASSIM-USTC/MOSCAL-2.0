#pragma once
#include "algebra.hpp"
#include "deom.hpp"
#include "linearalgebra.hpp"

void construct_Mapping_filter_p(DEOM_DATA *d, Trie *tree, int iado) {}

void construct_Mapping_filter_p(nNNmat &ddos, DEOM_DATA *d, Trie *tree,
                                int iado) {}

void rem_oprt_filter_p(nNNmat &ddos1, const nNNmat &total,
                       const DEOM_DATA *d, const int iado,
                       const char lcr) {}

void rem_oprt_1_filter_p(nNNmat &ddos1, const nNNmat &total,
                         const DEOM_DATA *d, const int iado,
                         const char lcr) {}

void rem_cal_filter_p(nNNmat &ddos1, const nNNmat &total,
                      const DEOM_DATA *d, const int iado) {}

void rem_cal_hei_filter_p(nNNmat &ddos1, const nNNmat &total,
                          const DEOM_DATA *d, const int iado) {}

void rem_cal_1_filter_p(nNNmat &ddos1, const nNNmat &total,
                        const DEOM_DATA *d, const int iado) {}

void rem_cal_1_hei_filter_p(nNNmat &ddos1, const nNNmat &total,
                            const DEOM_DATA *d, const int iado) {}

void syl_cal_ddos_filter_p(nNNmat &ddos1, const nNNmat &total,
                           const nNNmat &ddos2, const DEOM_DATA *d,
                           const DEOM_SYL *syl, const int iado) {}

void syl_cal_hei_filter_p(nNNmat &ddos1, const nNNmat &total,
                          const nNNmat &ddos2, const DEOM_DATA *d,
                          const DEOM_SYL *syl, const int iado) {}

void syl_cal_filter_p(nNNmat &ddos1, const nNNmat &total,
                      const DEOM_DATA *d, const DEOM_SYL *syl,
                      const int iado) {}

void syl_cal_1_filter_p(nNNmat &ddos1, const nNNmat &total,
                        const nNNmat &ddos2, const DEOM_DATA *d,
                        const DEOM_SYL *syl, const int iado) {}

void syl_cal_1_hei_filter_p(nNNmat &ddos1, const nNNmat &total,
                            const nNNmat &ddos2, const DEOM_DATA *d,
                            const DEOM_SYL *syl, const int iado) {}

void rem_cal_imag_filter_p(nNNmat &ddos1, const nNNmat &total,
                           const DEOM_DATA *d, const int iado) {}

void rem_cal_cfw_filter_p(nNNmat &ddos1, const nNNmat &total,
                          const DEOM_DATA *d, const int iado) {}

complex<double> caculate_trace(DEOM_DATA *d, Trie *tree) {
  MatrixXcd sum_Matrix = MatrixXcd::Zero(NSYS, NSYS);
  return sum_Matrix.trace();
}