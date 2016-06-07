#include<mpfr.h>
#include"Gfunction.h"
#include"theta.h"

void Gdelta::showLog(int base) {
  printf("Computing log h(n), it may be long ... \n");
  mpz_class num = Gprov.get_num();
  mpz_class den = Gprov.get_den();
  mpfr_set_default_prec(128);
  mpfr_t numer;
  mpfr_t denom;
  mpfr_t lognum;
  mpfr_t logden;
  mpfr_t res;
  mpfr_t thetapk;
  mpfr_inits(numer, denom, lognum, logden, res, thetapk, (mpfr_ptr) 0);
  mpfr_set_z(numer, num.get_mpz_t(), MPFR_RNDN);
  mpfr_set_z(denom, den.get_mpz_t(), MPFR_RNDN);
  mpfr_log(lognum, numer, MPFR_RNDN);
  mpfr_log(logden, denom, MPFR_RNDN);
  mpfr_sub(res, lognum, logden, MPFR_RNDN);
  //mpfr_printf("Glog: %.26Re\n\n",res);
  theta(thetapk, pk);
  //pfr_printf("Theta: %.26Re\n\n",thetapk);
  mpfr_add(res, res, thetapk, MPFR_RNDN);
  if (base==0)
    mpfr_printf("Log h(n): %.26Re\n\n",res);
  else {
    mpfr_t base_mpfr;
    mpfr_t base_log;
    mpfr_inits(base_mpfr, base_log, (mpfr_ptr) 0);
    mpfr_set_si(base_mpfr, base, MPFR_RNDN); 
    mpfr_log(base_log, base_mpfr, MPFR_RNDN);
    mpfr_div(res, res, base_log, MPFR_RNDN);
    mpfr_printf("Log h(n): %.26Re\n\n",res);
  }
}

