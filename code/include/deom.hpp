#pragma once
#include "headfile.hpp"
#include "pulse.hpp"

typedef struct {
  int nind;
  int lwsa;
  int lsad;
  MatrixXi twsa;
  MatrixXi tsad;
  ArrayXi dmax;
  ArrayXi amax;
} DEOM_SYM;

typedef struct {
  complex<double> Lambda_p;
  complex<double> Lambda_n;
  complex<double> tau;
  complex<double> alpha;
  complex<double> tf;
} CFW;

typedef struct {
  bool on;
  nNNmat sdip;
  nNNmat renormalize;
  VectorXcd bdip0;
  VectorXcd bdip1;
  VectorXcd bdip2;
  MatrixNcd pdip0;  // For single mode BSM-HQME
  MatrixNcd pdip1;  // For single mode BSM-HQME
  MatrixNcd pdip2;  // For single mode BSM-HQME
  nNNmat qmdta_r;
  nNNmat qmdta_l;
  nNNmat qmdtc_l;
  nNNmat qmdtc_r;
  nNNmat qmdt2a_l;
  nNNmat qmdt2a_r;
  nNNmat qmdt2b_l;
  nNNmat qmdt2b_r;
  nNNmat qmdt2c_l;
  nNNmat qmdt2c_r;
  nNNmat pdip;
  VectorXcd bdip;
} DIPOLE;

typedef struct {
  int flag_counter = 0;
  complex<double> trace = 0.0;
  complex<double> trace_1 = 0.0;
  bool flag = true;
} CONVERGED_FLAG;

typedef struct {
  int nmax;
  float ferr;
  int combmax;
  int nind;
  int lmax;
  int lmax_fp;
  int lmax_ma;
  int nmod;
  //
  MatrixNcd ham1;
  MatrixNcd hamt;
  double ham1_abs = 0.0;
  //
  MatrixNcd qmd1;  // For single mode BSM-HQME
  MatrixNcd qmdt;  // For single mode BSM-HQME
  nNNmat qmdta_l;
  nNNmat qmdta_r;
  nNNmat qmdtc_l;
  nNNmat qmdtc_r;
  VectorXd qmd1_abs;
  //
  double alp0;
  double alp1;
  double alp2;
  VectorXcd alp_list;
  int length_alp;
  // 
  nNNmat renormalize;
  nNNmat qmd2a;
  nNNmat qmd2b;
  nNNmat qmd2c;
  nNNmat qmdt2a_l;
  nNNmat qmdt2a_r;
  nNNmat qmdt2b_l;
  nNNmat qmdt2b_r;
  nNNmat qmdt2c_l;
  nNNmat qmdt2c_r;
  MatrixNb qmdt2a_equal_0;  // accelerate quadratic coupling
  MatrixNb qmdt2b_equal_0;  // accelerate quadratic coupling
  MatrixNb qmdt2c_equal_0;  // accelerate quadratic coupling
  //
  complex<double> omgb_BO;
  complex<double> one_over_2_times_mass;
  complex<double> coef0;  // constant for convenience
  complex<double> coef1;  // constant for convenience
  complex<double> coef2;  // constant for convenience
  complex<double> coef3;  // constant for convenience
  //
  DIPOLE dipole;
  DIPOLE dipole1;
  //
  VectorXcd expn;
  VectorXcd coef_abs;
  VectorXcd coef_lft;  // useless at the most cases
  VectorXcd coef_rht;  // useless at the most cases
  //
  MatrixNcd I;
  complex<double> i = complex<double>(0, 1);
  //
  int nddo;
  llmat comb_list;
  key_vec zerokey;
  key_vec emptykey;
  key_vec tempkey;
  DEOM_SYM sym;
  PULSE pulse;
  CFW cfw;
  //
  MatrixNcd V;
  MatrixNcd VI;
  VectorXcd S;
  //
  nNNmat ddos;
  nNNmat ddos1;
  nNNmat ddos2;
  nNNmat ddos3;
  nNNmat ddos4;
  nNNmat ddos5;
  nNNmat ddos6;
  nNNmat ddos7;
  nNNmat ddos8;
  nNNmat ddos9;
  nNNmat ddos10;
  array_key_vec keys;
  array_key_vec keys1;
  array_g_index g_index_list_pos;
  //
  vector<int> nsave = vector<int>(20);
  CONVERGED_FLAG if_converged;
  //
  Trie* tree;
  Trie* tree_hei;
  // some auxiliary json
  Json empty_list;
  Json empty_element;
} DEOM_DATA;

typedef struct {
  nNNmat ddos;
  int nddo;
} DEOM_BAK;

typedef struct {
  double dt;
  double ti;
  double tf;
  double dt2;
  double dt6;
  double ferr_adapt;
  char lcr;
  int nt_i;
  int nt_f;
  int nt;
  int step;
  double ferr;
  FILE* fpolar;
  FILE* frho;
  FILE* fcur;
  ArrayXd w_l;
  int w_len;
} CTRL;

typedef struct {
  int nind;
  int nddo;
  MatrixXi keys;
  MatrixNcd V;
  MatrixNcd VI;
  VectorXcd S;
  nNNmat U;
  nNNmat UI;
  nXvec R;
  double OMG;
  int lwsg;
  MatrixXi twsg;
  VectorXcd expn;
} DEOM_SYL;
