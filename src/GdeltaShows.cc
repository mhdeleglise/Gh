#include<gmpxx.h>
#include<gmp.h>
#include<mpfr.h>
#include"theta.h"
//#include<cmath>
//#include<vector>
//#include"numTheory.h"
//#include"utilities.h"
#include"Gfunction.h"
void prevprime(mpz_t rop, mpz_t x);
void Nk_compute(mpz_t Nk, long pk);

void Gdelta::show_value() {
  if (pk > 230561) {
    printf("The number of digits of h(n) is approximatively  %.1e\n", pk/::log((double) 10));///::log((double)10)).
    printf("I dont'think it is usefull to write this number.\n");
  }
 else
   {
     Nk_compute(Nk, pk);
     mpz_class numer=Gprov.get_num();
     mpz_class denom=Gprov.get_den();
     mpz_div(Nk, Nk, denom.get_mpz_t());
     mpz_mul(Nk, Nk, numer.get_mpz_t());
     gmp_printf("%.Zd\n",Nk);
   }
}


void Gdelta::show_factors() {
  int cnte=0;
  mpz_t q;
  mpz_init_set_si(q,pk);
  mpz_class numer=Gprov.get_num();
  mpz_class denom=Gprov.get_den();

  printf("[2-%ld]",pk);
  if (Gprov==1)
    {
      printf("\n");
      return;
    }
  
  printf(" ");
  while (numer > 1)
    {
      mpz_nextprime(q,q);
      while (! mpz_divisible_p(numer.get_mpz_t(), q)) {
	mpz_nextprime(q,q);
      }
      mpz_divexact(numer.get_mpz_t(), numer.get_mpz_t(), q);
      cnte+=1;
      if (cnte==1)
	if (numer==1) {
	  gmp_printf("x %.Zd ", q);
	  break;
	 }
	else {
	  gmp_printf("x %.Zd x ", q);
	}
      else {
	if (numer > 1)
	  gmp_printf("%.Zd x ", q);
	else
	  gmp_printf("%.Zd ", q);
      }
    }

  if (cnte==1) {
    gmp_printf(" / %.Zd\n",denom.get_mpz_t());
    return;
  }
  printf(" / ");  
  int i=0;
  while (i < cnte-1) {
    while (!mpz_divisible_p(denom.get_mpz_t(), q)) {
      prevprime(q,q);
    }
    i+=1;
    mpz_divexact(denom.get_mpz_t(), denom.get_mpz_t(), q);
    if (i < cnte-1)
      gmp_printf("%.Zd / ",q);
    else {
      gmp_printf("%.Zd / ",q);
      gmp_printf("%.Zd",denom.get_mpz_t());
    }
  }
  printf("\n");
}


void Gdelta::show_log() {
  cout << "Log   = " << log() << endl << endl;
}

void Gdelta::show_frac() {
  cout << "Gfrac = " << Gprov << endl << endl;;
}

void Gdelta::show_pk_Gp() {
  cout << pk << " " << Gprov.get_num() << " " << Gprov.get_den() << endl;
}


void Gdelta::showLog(int base) {
  //printf("Computing log h(n), it may be long .. a few seconds for n=10^20, \n");
  //printf("and the computaion time for h(10n) is twice the time of computation of h(n)\n");
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
  //  mpfr_printf("Glog: %.26Re\n\n",res);
  //cout << "Will compute theta pk= " << pk << endl;
  theta(thetapk, pk);
  //mpfr_printf("Theta: %.26Re\n\n",thetapk);
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

