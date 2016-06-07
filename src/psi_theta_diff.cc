#include<iostream>
#include<iomanip>
#include<mylib.h>
#include<mpfr.h>

mpfr_t _logp;
mpfr_t _pfloat;


void psi_theta_diff(mpfr_t res, long x) {
  mpfr_init_set_si(_logp, 0, MPFR_RNDN);
  mpfr_init_set_si(_pfloat, 0, MPFR_RNDN);
  mpfr_set_si(res, 0, MPFR_RNDN);
  int i=1;
  primes::init_primes(sqrt(x)+100);
  long p=primes::prime(i);
  long xp=x/p;
  while (p <= xp) {
    //cout << "p  " << p << "   xp " << xp << endl;
    mpfr_set_si(_pfloat, p, MPFR_RNDN);
    mpfr_log(_logp, _pfloat, MPFR_RNDN);
    while (p <= xp) {
      while (p <= xp) {
	mpfr_add(res, res, _logp, MPFR_RNDN);
	xp /= p;
      }
    }
    p = primes::prime(++i);
    xp = x/p;
  }
  mpfr_clear(_logp);
  mpfr_clear(_pfloat);
}

