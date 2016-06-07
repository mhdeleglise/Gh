#include <iostream>
#include <sstream>
#include <iomanip>
#include<mpfr.h>
#include<mylib.h>
#include "psi.h"

using namespace std;

void psi_theta_diff(mpfr_t res, long x);
  
void theta(mpfr_t res, long x) {
  mpfr_t delta;
  mpfr_init_set_si(delta,0, MPFR_RNDN);
  compute_psi(res, x);
  //mpfr_printf ("psi %.30RZf\n", res);
  psi_theta_diff(delta, x);
  //mpfr_printf ("delta= %.30RZf\n", delta);
  mpfr_sub(res, res, delta, MPFR_RNDN);
  mpfr_clear(delta);
}

//Computes theta(x) starting from the know value psi_value of psi(x)

void theta(mpfr_t res, long x, mpfr_t psi_value) {
  mpfr_t delta;
  mpfr_init_set_si(delta,0, MPFR_RNDN);
  psi_theta_diff(delta, x);
  //mpfr_printf ("delta= %.33RZf\n", delta);
  mpfr_set(res, psi_value, MPFR_RNDN);
  //mpfr_printf ("psi= %.33RZf\n", psi_value);
  //mpfr_printf ("res= %.33RZf\n", res);
  mpfr_sub(res, res, delta, MPFR_RNDN);
  mpfr_clear(delta);
}


