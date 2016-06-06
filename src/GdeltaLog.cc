#include<mpfr.h>
#include"Gfunction.h"

void Gdelta::showLog() {
  mpz_class num = Gprov.get_num();
  mpz_class den = Gprov.get_den();
  mpfr_set_default_prec(128);
  mpfr_t numer;
  mpfr_t denom;
  mpfr_t lognum;
  mpfr_t logden;
  mpfr_t res;
  mpfr_inits(numer, denom, lognum, logden, res, (mpfr_ptr) 0);
  mpfr_set_z(numer, num.get_mpz_t(), MPFR_RNDN);
  mpfr_set_z(denom, den.get_mpz_t(), MPFR_RNDN);
  mpfr_log(lognum, numer, MPFR_RNDN);
  mpfr_log(logden, denom, MPFR_RNDN);
  mpfr_sub(res, lognum, logden, MPFR_RNDN);
  mpfr_printf("Log: %.26Re\n\n",res);
}

