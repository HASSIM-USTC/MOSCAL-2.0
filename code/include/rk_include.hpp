#pragma once
#include "deom.hpp"

void RK1_filter(int i, int nsave1, int nsave2, nNNmat &ddos, nNNmat &ddos1,
                nNNmat &ddos3, double dt2) {
#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nsave1; i++) ddos3[i] = ddos[i] + ddos1[i] * dt2;

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave1; i < nsave2; i++) ddos3[i] = ddos1[i] * dt2;
}

void RK2_filter(int i, int nsave1, int nsave2, int nsave3, nNNmat &ddos,
                nNNmat &ddos1, nNNmat &ddos2, nNNmat &ddos3, double dt2) {
#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nsave1; i++) {
    ddos1[i] += ddos2[i] * 2.0;
    ddos3[i] = ddos[i] + ddos2[i] * dt2;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave1; i < nsave2; i++) {
    ddos1[i] += ddos2[i] * 2.0;
    ddos3[i] = ddos2[i] * dt2;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave2; i < nsave3; i++) {
    ddos1[i] = ddos2[i] * 2.0;
    ddos3[i] = ddos2[i] * dt2;
  }
}

void RK3_filter(int i, int nsave1, int nsave3, int nsave4, nNNmat &ddos,
                nNNmat &ddos1, nNNmat &ddos2, nNNmat &ddos3, double dt) {
#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nsave1; i++) {
    ddos1[i] += ddos2[i] * 2.0;
    ddos3[i] = ddos[i] + ddos2[i] * dt;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave1; i < nsave3; i++) {
    ddos1[i] += ddos2[i] * 2.0;
    ddos3[i] = ddos2[i] * dt;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave3; i < nsave4; i++) {
    ddos1[i] = ddos2[i] * 2.0;
    ddos3[i] = ddos2[i] * dt;
  }
}

void RK4_filter(int i, int nsave1, int nsave4, int nsave5, nNNmat &ddos,
                nNNmat &ddos1, nNNmat &ddos2, double dt6) {
#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nsave1; i++) {
    ddos1[i] += ddos2[i];
    ddos[i] += ddos1[i] * dt6;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave1; i < nsave4; i++) {
    ddos1[i] += ddos2[i];
    ddos[i] = ddos1[i] * dt6;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave4; i < nsave5; i++) {
    ddos1[i] = ddos2[i];
    ddos[i] = ddos1[i] * dt6;
  }
}

void RK1(int i, int nddo, nNNmat &ddos, nNNmat &ddos1, nNNmat &ddos3,
         double dt2) {
#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nddo; i++) ddos3[i] = ddos[i] + ddos1[i] * dt2;
}

void RK2(int i, int nddo, nNNmat &ddos, nNNmat &ddos1, nNNmat &ddos2,
         nNNmat &ddos3, double dt2) {
#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nddo; i++) {
    ddos1[i] += ddos2[i] * 2.0;
    ddos3[i] = ddos[i] + ddos2[i] * dt2;
  }
}

void RK3(int i, int nddo, nNNmat &ddos, nNNmat &ddos1, nNNmat &ddos2,
         nNNmat &ddos3, double dt) {
#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nddo; i++) {
    ddos1[i] += ddos2[i] * 2.0;
    ddos3[i] = ddos[i] + ddos2[i] * dt;
  }
}

void RK4(int i, int nddo, nNNmat &ddos, nNNmat &ddos1, nNNmat &ddos2,
         nNNmat &ddos3, double dt6) {
#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nddo; i++) {
    ddos1[i] += ddos2[i];
    ddos[i] += ddos1[i] * dt6;
  }
}

void RK32_1_filter(int i, int nsave1, int nsave2, nNNmat &ddos, nNNmat &ddos1,
                   nNNmat &ddos3, double dt) {
#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nsave1; i++) ddos3[i] = ddos[i] + ddos1[i] * dt / 2.0;

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave1; i < nsave2; i++) ddos3[i] = ddos1[i] * dt / 2.0;
}

void RK32_2_filter(int i, int nsave1, int nsave2, int nsave3, nNNmat &ddos,
                   nNNmat &ddos1, nNNmat &ddos2, nNNmat &ddos3, nNNmat &ddos4,
                   double dt) {
#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nsave1; i++) {
    ddos4[i] = ddos[i] + ddos1[i] * dt * 7.0 / 24.0 + ddos2[i] * dt / 4.0;
    ddos3[i] = ddos[i] + ddos2[i] * dt * 3.0 / 4.0;
    ddos1[i] = ddos[i] + ddos1[i] * dt * 2.0 / 9.0 + ddos2[i] * dt / 3.0;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave1; i < nsave2; i++) {
    ddos4[i] = ddos1[i] * dt * 7.0 / 24.0 + ddos2[i] * dt / 4.0;
    ddos3[i] = ddos2[i] * dt * 3.0 / 4.0;
    ddos1[i] = ddos1[i] * dt * 2.0 / 9.0 + ddos2[i] * dt / 3.0;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave2; i < nsave3; i++) {
    ddos4[i] = ddos2[i] * dt / 4.0;
    ddos3[i] = ddos2[i] * dt * 3.0 / 4.0;
    ddos1[i] = ddos2[i] * dt / 3.0;
  }
}

void RK32_3_filter(int i, int nsave1, int nsave3, int nsave4, nNNmat &ddos,
                   nNNmat &ddos1, nNNmat &ddos2, nNNmat &ddos3, nNNmat &ddos4,
                   double dt) {
#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nsave1; i++) {
    ddos3[i] = ddos1[i] + ddos2[i] * dt * 4.0 / 9.0;
    ddos4[i] += ddos2[i] * dt / 3.0;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave1; i < nsave3; i++) {
    ddos3[i] = ddos1[i] + ddos2[i] * dt * 4.0 / 9.0;
    ddos4[i] += ddos2[i] * dt / 3.0;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave3; i < nsave4; i++) {
    ddos3[i] = ddos2[i] * dt * 4.0 / 9.0;
    ddos4[i] = ddos2[i] * dt / 3.0;
  }
}

void RK32_4_filter(int i, int nsave4, int nsave5, nNNmat &ddos2, nNNmat &ddos3,
                   nNNmat &ddos4, double dt) {
#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nsave4; i++) {
    ddos4[i] += ddos2[i] * dt / 8.0 - ddos3[i];
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave4; i < nsave5; i++) {
    ddos4[i] = ddos2[i] * dt / 8.0;
  }
}

void adapt_rk(int i, int nsave4, int nsave5, nNNmat &ddos, nNNmat &ddos3,
              nNNmat &ddos4, double p, CTRL *c) {
  double ferr = 0;

#pragma omp parallel default(shared)
#pragma omp for reduction(max : ferr)
  for (i = 0; i < nsave5; i++) {
    for (int i_col = 0; i_col < NSYS; i_col++) {
      for (int i_row = 0; i_row < NSYS; i_row++) {
        ferr = max(ferr, abs((ddos4[i]).coeff(i_col, i_row)));
      }
    }
  }

  if (ferr > c->ferr_adapt) {
    c->dt = c->dt * 0.9 * pow(c->ferr_adapt / ferr, p);
  } else {
    c->ti += c->dt;
#pragma omp parallel default(shared)
#pragma omp for private(i) schedule(dynamic, 64)
    for (i = 0; i < nsave4; i++) ddos[i] = ddos3[i];
    c->dt = c->dt * 0.9 * pow(c->ferr_adapt / ferr, p);
  }
}

void RK54_1_filter(int i, vector<int> &nsave, nNNmat &ddos, nNNmat &ddos1,
                   nNNmat &ddos8, double dt) {
#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nsave[1]; i++) ddos8[i] = ddos[i] + ddos1[i] * dt / 5.0;

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave[1]; i < nsave[2]; i++) ddos8[i] = ddos1[i] * dt / 5.0;
}

void RK54_2_filter(int i, vector<int> &nsave, nNNmat &ddos, nNNmat &ddos1,
                   nNNmat &ddos2, nNNmat &ddos8, nNNmat &ddos9, nNNmat &ddos10,
                   double dt) {
#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nsave[1]; i++) {
    ddos8[i] =
        ddos[i] + ddos1[i] * dt * 3.0 / 40.0 + ddos2[i] * dt * 9.0 / 40.0;
    ddos9[i] = ddos[i] + ddos1[i] * dt * 35.0 / 384.0;
    ddos10[i] = ddos[i] + ddos1[i] * dt * 5179.0 / 57600.0;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave[1]; i < nsave[2]; i++) {
    ddos8[i] = ddos1[i] * dt * 3.0 / 40.0 + ddos2[i] * dt * 9.0 / 40.0;
    ddos9[i] = ddos1[i] * dt * 35.0 / 384.0;
    ddos10[i] = ddos1[i] * dt * 5179.0 / 57600.0;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave[2]; i < nsave[3]; i++) {
    ddos8[i] = ddos2[i] * dt * 9.0 / 40.0;
  }
}

void RK54_3_filter(int i, vector<int> &nsave, nNNmat &ddos, nNNmat &ddos1,
                   nNNmat &ddos2, nNNmat &ddos3, nNNmat &ddos8, nNNmat &ddos9,
                   nNNmat &ddos10, double dt) {
#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nsave[1]; i++) {
    ddos8[i] = ddos[i] + ddos1[i] * dt * 44.0 / 45.0 -
               ddos2[i] * dt * 56.0 / 15.0 + ddos3[i] * dt * 32.0 / 9.0;
    ddos9[i] = ddos9[i] + ddos3[i] * dt * 500.0 / 1113.0;
    ddos10[i] = ddos10[i] + ddos3[i] * dt * 7571.0 / 16695.0;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave[1]; i < nsave[2]; i++) {
    ddos8[i] = ddos1[i] * dt * 44.0 / 45.0 - ddos2[i] * dt * 56.0 / 15.0 +
               ddos3[i] * dt * 32.0 / 9.0;
    ddos9[i] = ddos9[i] + ddos3[i] * dt * 500.0 / 1113.0;
    ddos10[i] = ddos10[i] + ddos3[i] * dt * 7571.0 / 16695.0;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave[2]; i < nsave[3]; i++) {
    ddos8[i] = -ddos2[i] * dt * 56.0 / 15.0 + ddos3[i] * dt * 32.0 / 9.0;
    ddos9[i] = ddos3[i] * dt * 500.0 / 1113.0;
    ddos10[i] = ddos3[i] * dt * 7571.0 / 16695.0;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave[3]; i < nsave[4]; i++) {
    ddos8[i] = ddos3[i] * dt * 32.0 / 9.0;
    ddos9[i] = ddos3[i] * dt * 500.0 / 1113.0;
    ddos10[i] = ddos3[i] * dt * 7571.0 / 16695.0;
  }
}

void RK54_4_filter(int i, vector<int> &nsave, nNNmat &ddos, nNNmat &ddos1,
                   nNNmat &ddos2, nNNmat &ddos3, nNNmat &ddos4, nNNmat &ddos8,
                   nNNmat &ddos9, nNNmat &ddos10, double dt) {
#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nsave[1]; i++) {
    ddos8[i] = ddos[i] + ddos1[i] * dt * 19372.0 / 6561.0 -
               ddos2[i] * dt * 25360.0 / 2187.0 +
               ddos3[i] * dt * 64448.0 / 6561.0 - ddos4[i] * dt * 212.0 / 729.0;
    ddos9[i] = ddos9[i] + ddos4[i] * dt * 125.0 / 192.0;
    ddos10[i] = ddos10[i] + ddos4[i] * dt * 393.0 / 640.0;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave[1]; i < nsave[2]; i++) {
    ddos8[i] = ddos1[i] * dt * 19372.0 / 6561.0 -
               ddos2[i] * dt * 25360.0 / 2187.0 +
               ddos3[i] * dt * 64448.0 / 6561.0 - ddos4[i] * dt * 212.0 / 729.0;
    ddos9[i] = ddos9[i] + ddos4[i] * dt * 125.0 / 192.0;
    ddos10[i] = ddos10[i] + ddos4[i] * dt * 393.0 / 640.0;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave[2]; i < nsave[3]; i++) {
    ddos8[i] = -ddos2[i] * dt * 25360.0 / 2187.0 +
               ddos3[i] * dt * 64448.0 / 6561.0 - ddos4[i] * dt * 212.0 / 729.0;
    ddos9[i] = ddos9[i] + ddos4[i] * dt * 125.0 / 192.0;
    ddos10[i] = ddos10[i] + ddos4[i] * dt * 393.0 / 640.0;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave[3]; i < nsave[4]; i++) {
    ddos8[i] = ddos3[i] * dt * 64448.0 / 6561.0 - ddos4[i] * dt * 212.0 / 729.0;
    ddos9[i] = ddos9[i] + ddos4[i] * dt * 125.0 / 192.0;
    ddos10[i] = ddos10[i] + ddos4[i] * dt * 393.0 / 640.0;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave[4]; i < nsave[5]; i++) {
    ddos8[i] = -ddos4[i] * dt * 212.0 / 729.0;
    ddos9[i] = ddos4[i] * dt * 125.0 / 192.0;
    ddos10[i] = ddos4[i] * dt * 393.0 / 640.0;
  }
}

void RK54_5_filter(int i, vector<int> &nsave, nNNmat &ddos, nNNmat &ddos1,
                   nNNmat &ddos2, nNNmat &ddos3, nNNmat &ddos4, nNNmat &ddos5,
                   nNNmat &ddos8, nNNmat &ddos9, nNNmat &ddos10, double dt) {
#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nsave[1]; i++) {
    ddos8[i] = ddos[i] + ddos1[i] * dt * 9017.0 / 3168.0 -
               ddos2[i] * dt * 355.0 / 33.0 + ddos3[i] * dt * 46732.0 / 5247.0 +
               ddos4[i] * dt * 49.0 / 176.0 - ddos5[i] * dt * 5103.0 / 18656.0;
    ddos9[i] = ddos9[i] - ddos5[i] * dt * 2187.0 / 6784.0;
    ddos10[i] = ddos10[i] - ddos5[i] * dt * 92097.0 / 339200.0;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave[1]; i < nsave[2]; i++) {
    ddos8[i] = ddos1[i] * dt * 9017.0 / 3168.0 - ddos2[i] * dt * 355.0 / 33.0 +
               ddos3[i] * dt * 46732.0 / 5247.0 + ddos4[i] * dt * 49.0 / 176.0 -
               ddos5[i] * dt * 5103.0 / 18656.0;
    ddos9[i] = ddos9[i] - ddos5[i] * dt * 2187.0 / 6784.0;
    ddos10[i] = ddos10[i] - ddos5[i] * dt * 92097.0 / 339200.0;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave[2]; i < nsave[3]; i++) {
    ddos8[i] = -ddos2[i] * dt * 355.0 / 33.0 +
               ddos3[i] * dt * 46732.0 / 5247.0 + ddos4[i] * dt * 49.0 / 176.0 -
               ddos5[i] * dt * 5103.0 / 18656.0;
    ddos9[i] = ddos9[i] - ddos5[i] * dt * 2187.0 / 6784.0;
    ddos10[i] = ddos10[i] - ddos5[i] * dt * 92097.0 / 339200.0;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave[3]; i < nsave[4]; i++) {
    ddos8[i] = ddos3[i] * dt * 46732.0 / 5247.0 + ddos4[i] * dt * 49.0 / 176.0 -
               ddos5[i] * dt * 5103.0 / 18656.0;
    ddos9[i] = ddos9[i] - ddos5[i] * dt * 2187.0 / 6784.0;
    ddos10[i] = ddos10[i] - ddos5[i] * dt * 92097.0 / 339200.0;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave[4]; i < nsave[5]; i++) {
    ddos8[i] = ddos4[i] * dt * 49.0 / 176.0 - ddos5[i] * dt * 5103.0 / 18656.0;
    ddos9[i] = ddos9[i] - ddos5[i] * dt * 2187.0 / 6784.0;
    ddos10[i] = ddos10[i] - ddos5[i] * dt * 92097.0 / 339200.0;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave[5]; i < nsave[6]; i++) {
    ddos8[i] = -ddos5[i] * dt * 5103.0 / 18656.0;
    ddos9[i] = -ddos5[i] * dt * 2187.0 / 6784.0;
    ddos10[i] = -ddos5[i] * dt * 92097.0 / 339200.0;
  }
}

void RK54_6_filter(int i, vector<int> &nsave, nNNmat &ddos6, nNNmat &ddos9,
                   nNNmat &ddos10, double dt) {
#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nsave[6]; i++) {
    ddos9[i] = ddos9[i] + ddos6[i] * dt * 11.0 / 84.0;
    ddos10[i] = ddos10[i] + ddos6[i] * dt * 187.0 / 2100.0;
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave[6]; i < nsave[7]; i++) {
    ddos9[i] = ddos6[i] * dt * 11.0 / 84.0;
    ddos10[i] = ddos6[i] * dt * 187.0 / 2100.0;
  }
}

void RK54_7_filter(int i, vector<int> &nsave, nNNmat &ddos7, nNNmat &ddos9,
                   nNNmat &ddos10, double dt) {
#pragma omp for private(i) schedule(dynamic, 64)
  for (i = 0; i < nsave[7]; i++) {
    ddos10[i] = ddos10[i] + ddos7[i] * dt * 1.0 / 40.0 - ddos9[i];
  }

#pragma omp for private(i) schedule(dynamic, 64)
  for (i = nsave[7]; i < nsave[8]; i++) {
    ddos10[i] = ddos7[i] * dt * 1.0 / 40.0;
  }
}
