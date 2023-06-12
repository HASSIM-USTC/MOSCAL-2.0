#pragma once
#include "deom.hpp"

void lambda_forward(DEOM_DATA *d, double t) {
#pragma omp single
  {
    d->cfw.Lambda_p = (1.0 - exp(-d->cfw.alpha * t)) / (1.0 - exp(-d->cfw.alpha * d->cfw.tf)) + d->cfw.tau / 2.0 * (d->cfw.alpha * exp(-d->cfw.alpha * t)) / (1.0 - exp(-d->cfw.alpha * d->cfw.tf));
    d->cfw.Lambda_n = (1.0 - exp(-d->cfw.alpha * t)) / (1.0 - exp(-d->cfw.alpha * d->cfw.tf)) - d->cfw.tau / 2.0 * (d->cfw.alpha * exp(-d->cfw.alpha * t)) / (1.0 - exp(-d->cfw.alpha * d->cfw.tf));
  }
}

void lambda_backward(DEOM_DATA *d, double t) {
#pragma omp single
  {
    d->cfw.Lambda_p = (exp(d->cfw.alpha * d->cfw.tf) - exp(d->cfw.alpha * t)) / (exp(d->cfw.alpha * d->cfw.tf) - 1.0) - d->cfw.tau / 2.0 * (d->cfw.alpha * exp(d->cfw.alpha * t)) / (exp(d->cfw.alpha * d->cfw.tf) - 1.0);
    d->cfw.Lambda_n = (exp(d->cfw.alpha * d->cfw.tf) - exp(d->cfw.alpha * t)) / (exp(d->cfw.alpha * d->cfw.tf) - 1.0) + d->cfw.tau / 2.0 * (d->cfw.alpha * exp(d->cfw.alpha * t)) / (exp(d->cfw.alpha * d->cfw.tf) - 1.0);
  }
}