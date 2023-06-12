#pragma once
#include "deom.hpp"

void print(FILE *fpolar, complex<double> trace, double t) {
  printf("%16.6e%20.10e%20.10e\n", t, real(trace), imag(trace));
  fprintf(fpolar, "%16.6e%20.10e%20.10e\n", t, real(trace), imag(trace));
  fflush(fpolar);
  fflush(stdout);
}

void print(complex<double> trace, double t) {
  printf("%16.6e%20.10e%20.10e\n", t, real(trace), imag(trace));
  fflush(stdout);
}

void print(FILE *frho, FILE *fcur, int ii, DEOM_DATA *d, const CTRL *c) {
  print_is_hermite(d->ddos[0], 10 * d->ferr);
  printf("%20.16e\t", c->ti + ii * c->dt);
  for (int i = 0; i < NSYS; i++)
    printf("%20.16e\t", d->ddos[0].coeff(i, i).real());
  printf("\n");

  fprintf(frho, "%20.16e\t", c->ti + ii * c->dt);
  for (int i = 0; i < NSYS; i++)
    for (int j = 0; j < NSYS; j++)
      fprintf(frho, "%20.16e\t%20.16e\t", d->ddos[0].coeff(i, j).real(),
              d->ddos[0].coeff(i, j).imag());
  fprintf(frho, "\n");

  fprintf(fcur, "%20.16e\t", c->ti);
  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    int pos = generate_key_plus(d->zerokey, d->emptykey, i_nind, d) - 1;
    fprintf(fcur, "%20.16e\t%20.16e\t",
            (d->i * MatrixXcd(d->qmdta_l[i_nind] * d->ddos[pos]).trace())
                .real(),
            (d->i * MatrixXcd(d->qmdta_l[i_nind] * d->ddos[pos]).trace())
                .imag());
  }
  fprintf(fcur, "\n");
  fflush(frho);
  fflush(stdout);
}

void print_md(FILE *frho, FILE *fcur, int ii, DEOM_DATA *d,
                     const CTRL *c) {
  printf("%20.16e\t", c->ti + (ii + 1) * c->dt);
  fprintf(frho, "%20.16e\t", c->ti + (ii + 1) * c->dt);

  for (int i_nind = 0; i_nind < 2; i_nind++) {
    int pos = generate_key_plus(d->zerokey, d->emptykey, i_nind, d) - 1;
    for (int i = 0; i < NSYS; i++) {
      printf("%20.16e\t", d->ddos[pos].coeff(i, i).real());
      fprintf(frho, "%20.16e\t%20.16e\t", d->ddos[pos].coeff(i, i).real(),
              d->ddos[pos].coeff(i, i).imag());
    }
  }
  printf("\n");
  fprintf(frho, "\n");
  fflush(stdout);
  fflush(frho);
}

void print_md_filter(FILE *frho, FILE *fcur, int ii, DEOM_DATA *d,
                     const CTRL *c) {
  printf("%20.16e\t", c->ti + (ii + 1) * c->dt);
  fprintf(frho, "%20.16e\t", c->ti + (ii + 1) * c->dt);

  for (int i_nind = 0; i_nind < 2; i_nind++) {
    generate_key_plus(d->zerokey, d->emptykey, i_nind, d);
    int pos = d->tree->find(gen_hash_value(d->emptykey, d));
    if (pos != -1) {
      for (int i = 0; i < NSYS; i++) {
        printf("%20.16e\t", d->ddos[pos].coeff(i, i).real());
        fprintf(frho, "%20.16e\t%20.16e\t",
                d->ddos[pos].coeff(i, i).real(),
                d->ddos[pos].coeff(i, i).imag());
      }
    } else {
      for(int i = 0; i < NSYS; i++) {
        fprintf(frho, "%20.16e\t%20.16e\t", 0.0, 0.0);
        printf("%20.16e\t", 0.0);
      }
    }
  }
  printf("\n");
  fprintf(frho, "\n");
  fflush(stdout);
  fflush(frho);
}

void print_filter(FILE *frho, FILE *fcur, int ii, DEOM_DATA *d,
                  const CTRL *c) {
  int pos = d->tree->find(1);

  print_is_hermite(d->ddos[pos], 10 * d->ferr);
  printf("%20.16e\t", c->ti + (ii + 1) * c->dt);
  for (int i = 0; i < NSYS; i++)
    printf("%20.16e\t", d->ddos[pos].coeff(i, i).real());
  printf("\n");

  fprintf(frho, "%20.16e\t", c->ti + (ii + 1) * c->dt);
  for (int i = 0; i < NSYS; i++)
    for (int j = 0; j < NSYS; j++)
      fprintf(frho, "%20.16e\t%20.16e\t", d->ddos[pos].coeff(i, j).real(),
              d->ddos[pos].coeff(i, j).imag());
  fprintf(frho, "\n");

  fprintf(fcur, "%20.16e\t", c->ti + (ii + 1) * c->dt);
  for (int i_nind = 0; i_nind < d->nind; i_nind++) {
    generate_key_plus(d->zerokey, d->emptykey, i_nind, d);
    pos = d->tree->find(gen_hash_value(d->emptykey, d));
    if (pos != -1) {
      fprintf(fcur, "%20.16e\t%20.16e\t",
              (d->i * MatrixXcd(d->qmdta_l[i_nind] * d->ddos[pos]).trace())
                  .real(),
              (d->i * MatrixXcd(d->qmdta_l[i_nind] * d->ddos[pos]).trace())
                  .imag());
    } else {
      fprintf(fcur, "%20.16e\t%20.16e\t", 0.0, 0.0);
    }
  }
  fprintf(fcur, "\n");
  fflush(stdout);
  fflush(frho);
}

void print_adapt_rk(FILE *frho, DEOM_DATA *d, CTRL *c) {
  int pos = d->tree->find(1);
  print_is_hermite(d->ddos[pos], 10 * d->ferr);

  printf("%20.16e\t", c->ti);
  for (int i = 0; i < NSYS; i++)
    printf("%20.16e\t", d->ddos[pos].coeff(i, i).real());
  printf("\n");

  fprintf(frho, "%20.16e\t", c->ti);
  for (int i = 0; i < NSYS; i++)
    for (int j = 0; j < NSYS; j++)
      fprintf(frho, "%20.16e\t%20.16e\t", d->ddos[pos].coeff(i, j).real(),
              d->ddos[pos].coeff(i, j).imag());
  fprintf(frho, "\n");
  fflush(stdout);
  fflush(frho);
}