#pragma once
#include "headfile.hpp"

typedef struct {
  bool on;
  double epsilon0;
  double delta0;
  double g0;
} PULSE;

MatrixNcd pulse_et_ham1(double t, PULSE p) {
  MatrixXcd _et;
  MatrixNcd _et_out;
  _et.setZero();
  double t1 = 3.14159265358979323846 / (4 * p.epsilon0);
  double t2 = t1 + 3.14159265358979323846 / (2 * p.g0);
  double t3 = t2 + 3.14159265358979323846 / (4 * p.delta0);
  double t4 = t3 + 3.14159265358979323846 / (2 * p.g0);
  double t5 = t4 + 3.14159265358979323846 / (2 * p.delta0);
  double t6 = t5 + 3.14159265358979323846 / (4 * p.epsilon0);
  double t7 = t6 + 3.14159265358979323846 / (4 * p.delta0);
  if (p.on) {
    if (t <= t1 && t >= 0) {
      _et(0, 0) = p.epsilon0 * 2;
      _et(3, 3) = -p.epsilon0 * 2;
    }
    if (t <= t2 && t > t1) _et(1, 2) = _et(2, 1) = -p.g0;
    if (t <= t3 && t > t2) {
      _et(0, 2) = _et(2, 0) = p.delta0;
      _et(1, 3) = _et(3, 1) = p.delta0;
    }
    if (t <= t4 && t > t3) _et(1, 2) = _et(2, 1) = p.g0;
    if (t <= t5 && t > t4) {
      _et(0, 1) = _et(1, 0) = p.delta0;
      _et(2, 3) = _et(3, 2) = p.delta0;
    }
    if (t <= t6 && t > t5) {
      _et(0, 0) = _et(2, 2) = p.epsilon0;
      _et(1, 1) = _et(3, 3) = -p.epsilon0;
    }
    if (t <= t7 && t > t6) {
      _et(0, 1) = _et(1, 0) = -p.delta0;
      _et(2, 3) = _et(3, 2) = -p.delta0;
    }
  }
  allocate_pulse(_et_out, _et);
  return _et_out;
}

void print(const PULSE pulse) {
  printf("epsilon0\t%e\n", pulse.epsilon0);
  printf("delta0\t%e\n", pulse.delta0);
  printf("g0\t%e\n", pulse.g0);
}

void Init_pulse(PULSE &pulse, const Json &json) {
  // ## init pulse##
  pulse.on = true;
  pulse.epsilon0 = json["epsilon0"];
  pulse.delta0 = json["delta0"];
  pulse.g0 = json["g0"];
  // ## init pulse##
}