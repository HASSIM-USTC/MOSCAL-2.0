#include "equilibrium.hpp"
#include "init.hpp"
// #include "read_write.hpp"
#include "thermo.hpp"

void init_deom(const Json &json, int argc, char *argv[]) {
  DEOM_DATA *d = new DEOM_DATA();
  CTRL *c = new CTRL();
  Init_sys(d, json);

  Init_ctrl(c, json["equilibrium"]);
  c->frho = fopen("prop-rho-eq1.dat", "w");
  c->fcur = fopen("curr.dat", "w");
  Init_sys_filter(d, json);
  Init_aux_dt(d, json);

#if defined(STD)
  omp_set_num_threads(1);
#endif  // STD

  if (json.value("cfw", false)) {
    Init_cfw(d, json["cfw-data"]);
    if (json["cfw-data"]["forward"]) {
      equilibrium_cfw_forward_filter(d, c);
    } else {
      fclose(c->frho);
      c->frho = fopen("prop-rho-eq.dat", "w");
      Init_ctrl(c, json["equilibrium-dt"]);
      equilibrium_filter(d, c);
      fclose(c->frho);
      c->frho = fopen("prop-rho-eq1.dat", "w");
      Init_ctrl(c, json["equilibrium"]);
      equilibrium_cfw_backward_filter(d, c);
    }
  }

  if (json.value("imag", false)) {
    equilibrium_imag_filter(d, c);
  }

  fclose(c->frho);
  fclose(c->fcur);
}

int main(int argc, char *argv[]) {
  ifstream jsonFile("input.json");
  Json json;
  jsonFile >> json;
  init_deom(json, argc, argv);
}
