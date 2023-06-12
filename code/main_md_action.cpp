#include "1d_corr.hpp"
#include "equilibrium.hpp"
#include "init.hpp"

void init_deom(const Json &json, int argc, char *argv[]) {
  DEOM_DATA *d = new DEOM_DATA();
  CTRL *c = new CTRL();
  Init_sys(d, json);

  if (json["filter"]) {
#if defined(STD)
    omp_set_num_threads(1);
#endif  // STD
    Init_ctrl(c, json["equilibrium"]);
    c->frho = fopen("prop-rho-eq.dat", "w");
    c->fcur = fopen("curr.dat", "w");
    Init_sys_filter(d, json);
    if (json["equilibrium"]["dt-method"]) {
      Init_aux_dt(d, json);
      equilibrium_md_filter(d, c);
    }
  } else {
    Init_ctrl(c, json["equilibrium"]);
    c->frho = fopen("prop-rho-eq.dat", "w");
    c->fcur = fopen("curr.dat", "w");
    if (json["equilibrium"]["dt-method"]) {
      Init_aux_dt(d, json);
      equilibrium_md(d, c);
    }
    // TODO destruction.
  }
}

int main(int argc, char *argv[]) {
  ifstream jsonFile("input.json");
  Json json;
  jsonFile >> json;
  init_deom(json, argc, argv);
}
