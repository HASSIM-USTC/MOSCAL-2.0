#include "1d_corr.hpp"
#include "equilibrium.hpp"
#include "init.hpp"

void calculator_deom(const Json &json, DEOM_DATA *d, CTRL *c,
                     DEOM_SYL *syl) {
  if (json["if-time"].get<bool>()) {
    Init_aux_dt(d, json);
    Init_ctrl(c, json["time"]);
    bool if_noise = json["time"].value("noise", true);
    if (if_noise)
      od_corr_noise(d, c);
    else
      od_corr(d, c);
  } else {
    Init_aux_spe_sc2(d, json);
    Init_syl_omg(syl, d, json["frequency"]);
    Init_ctrl(c, json["frequency"]);
    bool if_noise_w = json["frequency"].value("noise", true);
    if (if_noise_w)
      od_corr_sc2_noise(d, c, syl);
    else
      od_corr_sc2(d, c, syl);
  }
}

void calculator_deom_filter(const Json &json, DEOM_DATA *d, CTRL *c,
                            DEOM_SYL *syl) {
  if (json["if-time"].get<bool>()) {
    Init_aux_dt(d, json);
    Init_ctrl(c, json["time"]);
    bool if_noise = json["time"].value("noise", true);
    if (json["time"].value("Hei", false)) {
      if (if_noise)
        od_corr_noise_hei_filter(d, c);
      else
        od_corr_hei_filter(d, c);
    } else {
      if (if_noise)
        od_corr_noise_filter(d, c);
      else
        od_corr_filter(d, c);
    }
  } else {
    Init_aux_dt(d, json);
    Init_syl_omg(syl, d, json["frequency"]);
    Init_ctrl(c, json["frequency"]);
    bool if_noise_w = json["frequency"].value("noise", true);
    if (json["frequency"].value("Hei", false)) {
      if (if_noise_w)
        od_corr_Hei_sc2_noise_filter(d, c, syl);
      else
        od_corr_Hei_sc2_filter(d, c, syl);
    } else {
      if (if_noise_w)
        od_corr_sc2_noise_filter(d, c, syl);
      else
        od_corr_sc2_filter(d, c, syl);
    }
  }
}

template <typename FUN>
void string_script(const Json &json, DEOM_DATA *d, DEOM_BAK *b, CTRL *c,
                   DEOM_SYL *syl, string spectrum, string data,
                   FUN calculator_deom) {
  if (json.value(spectrum, false)) {
    Init_dip(&d->dipole, d, json[data]["dipole"]);
    Init_dip(&d->dipole1, d, json[data]["dipole1"]);
    c->fpolar = fopen(json[data]["file"].get<string>().c_str(), "w");
    calculator_deom(json[data], d, c, syl);
    fclose(c->fpolar);
  }

  if (json.value("backup", false)) {
    Init_bak(d, b);
  } else {
    exit(1);
  }
}

void init_deom(const Json &json, int argc, char *argv[]) {
  DEOM_DATA *d = new DEOM_DATA();
  CTRL *c = new CTRL();
  DEOM_SYL *syl = new DEOM_SYL();
  Init_sys(d, json);

  if (json["filter"].get<bool>()) {
#if defined(STD)
    omp_set_num_threads(1);
#endif  // STD

    Init_ctrl(c, json["equilibrium"]);
    c->frho = fopen("prop-rho-eq.dat", "w");
    c->fcur = fopen("curr.dat", "w");
    Init_sys_filter(d, json);
    if (json["equilibrium"].value("sc2", false)) {
      Init_aux_sc2(d, json);
      Init_syl(syl, d, json["syl"]);
      Init_syl_omg(syl, d, json["equilibrium"]);
      equilibrium_sc2_filter(d, c, syl);
    } else if (json["equilibrium"].value("dt-method", false)) {
      if (json["equilibrium"].value("dt-method-RK32", false)) {
        Init_aux_dt_32(d, json);
        equilibrium_32_filter(d, c);
      } else if (json["equilibrium"].value("dt-method-RK54", false)) {
        Init_aux_dt_54(d, json);
        equilibrium_54_filter(d, c);
      } else {
        Init_aux_dt(d, json);
        equilibrium_filter(d, c);
      }
    }

    fclose(c->frho);
    fclose(c->fcur);

    DEOM_BAK *b = new DEOM_BAK;
    if (json.value("backup", false)) {
      Init_bak(b, d);
    }

    string_script(json, d, b, c, syl, "spectrum", "spectrum-data",
                  calculator_deom_filter);
  } else {
    c->frho = fopen("prop-rho-eq.dat", "w");
    c->fcur = fopen("curr.dat", "w");
    Init_syl(syl, d, json["syl"]);
    Init_ctrl(c, json["equilibrium"]);
    Init_syl_omg(syl, d, json["equilibrium"]);
    if (json["equilibrium"]["sc2"].get<bool>()) {
      Init_aux_sc2(d, json);
      equilibrium_sc2(d, c, syl);
    } else if (json["equilibrium"]["dt-method"].get<bool>()) {
      Init_aux_dt(d, json);
      equilibrium(d, c);
    }

    fclose(c->frho);
    fclose(c->fcur);

    DEOM_BAK *b = new DEOM_BAK;
    if (json.value("backup", false)) {
      Init_bak(b, d);
    }

    string_script(json, d, b, c, syl, "spectrum", "spectrum-data",
                  calculator_deom);
    string_script(json, d, b, c, syl, "spectrum1", "spectrum-data1",
                  calculator_deom);
    string_script(json, d, b, c, syl, "spectrum2", "spectrum-data2",
                  calculator_deom);
    string_script(json, d, b, c, syl, "spectrum3", "spectrum-data3",
                  calculator_deom);
    string_script(json, d, b, c, syl, "spectrum4", "spectrum-data4",
                  calculator_deom);
    string_script(json, d, b, c, syl, "spectrum5", "spectrum-data5",
                  calculator_deom);
    string_script(json, d, b, c, syl, "spectrum6", "spectrum-data6",
                  calculator_deom);
    string_script(json, d, b, c, syl, "spectrum7", "spectrum-data7",
                  calculator_deom);
  }
  // TODO destruction.
}

int main(int argc, char *argv[]) {
  ifstream jsonFile("input.json");
  Json json;
  jsonFile >> json;
  init_deom(json, argc, argv);
}
