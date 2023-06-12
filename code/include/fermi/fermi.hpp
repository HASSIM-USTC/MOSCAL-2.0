#pragma once
#include "mr_fermi_gen.hpp"
#include "mr_fermi_p.hpp"

void gen_g_index(const DEOM_DATA *d, key_vec &key, ArrayXi &g_index_list,
                 int iado) {
  g_index_list(0) = iado + 1;
  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    if (key(i_nind)) {
      g_index_list(i_nind + 1) = generate_hash_value_minus(key, i_nind, d);
    } else if (tier(key, d) < d->lmax)
      g_index_list(i_nind + 1) = generate_hash_value_plus(key, i_nind, d);
    else
      g_index_list(i_nind + 1) = 0;
  }
}

void gen_g_index_filter(const DEOM_DATA *d, key_vec &key,
                        ArrayXi &g_index_list, int iado) {
  g_index_list = d->g_index_list_pos[iado];
}

void rem_oprt_filter_p(nNNmat &ddos1, const nNNmat &total,
                       const DEOM_DATA *d, const int iado,
                       const char lcr) {
  rem_oprt_gen(ddos1, total, d, iado, lcr, gen_g_index_filter);
}

void rem_oprt_p(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
                const int iado, const char lcr) {
  rem_oprt_gen(ddos1, total, d, iado, lcr, gen_g_index);
}

void rem_cal_filter_p(nNNmat &ddos1, const nNNmat &total,
                      const DEOM_DATA *d, const int iado) {
  rem_cal_gen(ddos1, total, d, iado, gen_g_index_filter);
}

void rem_cal_p(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
               const int iado) {
  rem_cal_gen(ddos1, total, d, iado, gen_g_index);
}

void rem_cal_hei_filter_p(nNNmat &ddos1, const nNNmat &total,
                          const DEOM_DATA *d, const int iado) {
  rem_cal_hei_gen(ddos1, total, d, iado, gen_g_index_filter);
}

void rem_cal_hei_p(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
                   const int iado) {
  rem_cal_hei_gen(ddos1, total, d, iado, gen_g_index);
}

void rem_cal_1_filter_p(nNNmat &ddos1, const nNNmat &total,
                        const DEOM_DATA *d, const int iado) {
  rem_cal_1_gen(ddos1, total, d, iado, gen_g_index_filter);
}

void rem_cal_1_p(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
                 const int iado) {
  rem_cal_1_gen(ddos1, total, d, iado, gen_g_index);
}

void rem_cal_1_hei_filter_p(nNNmat &ddos1, const nNNmat &total,
                            const DEOM_DATA *d, const int iado) {
  rem_cal_1_hei_gen(ddos1, total, d, iado, gen_g_index_filter);
}

void rem_cal_1_hei_p(nNNmat &ddos1, const nNNmat &total,
                     const DEOM_DATA *d, const int iado) {
  rem_cal_1_hei_gen(ddos1, total, d, iado, gen_g_index);
}

void syl_cal_ddos_filter_p(nNNmat &ddos1, const nNNmat &total,
                           const nNNmat &ddos2, const DEOM_DATA *d,
                           const DEOM_SYL *syl, const int iado) {
  syl_cal_ddos_gen(ddos1, total, ddos2, d, syl, iado, gen_g_index_filter);
}

void syl_cal_ddos_p(nNNmat &ddos1, const nNNmat &total,
                    const nNNmat &ddos2, const DEOM_DATA *d,
                    const DEOM_SYL *syl, const int iado) {
  syl_cal_ddos_gen(ddos1, total, ddos2, d, syl, iado, gen_g_index);
}

void syl_cal_filter_p(nNNmat &ddos1, const nNNmat &total,
                      const DEOM_DATA *d, const DEOM_SYL *syl,
                      const int iado) {
  syl_cal_gen(ddos1, total, d, syl, iado, gen_g_index_filter);
}

void syl_cal_p(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
               const DEOM_SYL *syl, const int iado) {
  syl_cal_gen(ddos1, total, d, syl, iado, gen_g_index);
}

void syl_cal_hei_filter_p(nNNmat &ddos1, const nNNmat &total,
                          const nNNmat &ddos2, const DEOM_DATA *d,
                          const DEOM_SYL *syl, const int iado) {
  syl_cal_Hei_gen(ddos1, total, ddos2, d, syl, iado, gen_g_index_filter);
}

void syl_cal_Hei_p(nNNmat &ddos1, const nNNmat &total, const nNNmat &ddos2,
                   const DEOM_DATA *d, const DEOM_SYL *syl,
                   const int iado) {
  syl_cal_Hei_gen(ddos1, total, ddos2, d, syl, iado, gen_g_index);
}

void syl_cal_1_filter_p(nNNmat &ddos1, const nNNmat &total,
                        const nNNmat &ddos2, const DEOM_DATA *d,
                        const DEOM_SYL *syl, const int iado) {
  syl_cal_1_gen(ddos1, total, ddos2, d, syl, iado, gen_g_index_filter);
}

void syl_cal_1_p(nNNmat &ddos1, const nNNmat &total, const nNNmat &ddos2,
                 const DEOM_DATA *d, const DEOM_SYL *syl, const int iado) {
  syl_cal_1_gen(ddos1, total, ddos2, d, syl, iado, gen_g_index);
}

void syl_cal_1_hei_filter_p(nNNmat &ddos1, const nNNmat &total,
                            const nNNmat &ddos2, const DEOM_DATA *d,
                            const DEOM_SYL *syl, const int iado) {
  syl_cal_1_Hei_gen(ddos1, total, ddos2, d, syl, iado, gen_g_index_filter);
}

void syl_cal_1_Hei_p(nNNmat &ddos1, const nNNmat &total,
                     const nNNmat &ddos2, const DEOM_DATA *d,
                     const DEOM_SYL *syl, const int iado) {
  syl_cal_1_Hei_gen(ddos1, total, ddos2, d, syl, iado, gen_g_index);
}

void rem_cal_imag_filter_p(nNNmat &ddos1, const nNNmat &total,
                           const DEOM_DATA *d, const int iado) {
  rem_cal_imag_gen(ddos1, total, d, iado, gen_g_index_filter);
}

void rem_cal_imag_p(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
                    const int iado) {
  rem_cal_imag_gen(ddos1, total, d, iado, gen_g_index);
}

void rem_cal_cfw_filter_p(nNNmat &ddos1, const nNNmat &total,
                          const DEOM_DATA *d, const int iado) {
  rem_cal_cfw_gen(ddos1, total, d, iado, gen_g_index_filter);
}

void rem_cal_cfw_p(nNNmat &ddos1, const nNNmat &total, const DEOM_DATA *d,
                   const int iado) {
  rem_cal_cfw_gen(ddos1, total, d, iado, gen_g_index);
}