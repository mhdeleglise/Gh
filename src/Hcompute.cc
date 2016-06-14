#include<gmpxx.h>
#include<cmath>
#include<vector>
#include<unistd.h>
#include"utilities.h"
#include"Gfunction.h"
#include"theta.h"
#include<mpfr.h>
#include<gmp.h>
#include"Nk.h"
void invli(mpfr_t res, mpfr_t x, double precision);

class Hn {
public:
  Gdelta gdelta;
  mpz_t n;
  mpfr_t n_mpfr;
  mpfr_t logn;
  mpfr_t thetak;
  mpfr_t logh;
  mpfr_t num_bn;
  mpfr_t den_bn;
  mpfr_t bn;
  mpfr_t invli2n;

  Hn(mpz_t n, long pk, long deltak);
  void show_Bvalue() {};
  void compute_Log();
  void compute_bn();
};

Hn::Hn(mpz_t nn, long pk,long m){
  gdelta.init(pk,m,500,0);
  mpz_init(n);
  mpz_set(n,nn);
  mpfr_init_set_z(n_mpfr,n,MPFR_RNDN);
  mpfr_inits2(256, logn, thetak, logh, num_bn, den_bn, bn, invli2n, (mpfr_ptr) 0);
  mpfr_set_si(logn,0,MPFR_RNDN);
  mpfr_set_si(thetak,0,MPFR_RNDN);
  mpfr_set_si(logh,0,MPFR_RNDN);
  mpfr_set_si(num_bn,0,MPFR_RNDN);
  mpfr_set_si(den_bn,0,MPFR_RNDN);
  mpfr_set_si(bn,0,MPFR_RNDN);
  mpfr_set_si(invli2n,0,MPFR_RNDN);
}

void Hn::compute_Log() {
  theta(thetak,gdelta.pk);
  gdelta.compute_Glog();
  mpfr_add(logh, thetak, gdelta.Glog, MPFR_RNDN);
}

void Hn::compute_bn() {
  mpfr_t invli2n,x2;
  mpfr_inits2(256,x2, invli2n, (mpfr_ptr) 0);
  mpfr_set_si(x2,0,MPFR_RNDN);

  invli(x2,n_mpfr,1e-12);
  mpfr_sqrt(invli2n,x2,MPFR_RNDN);
  compute_Log();
  mpfr_sub(num_bn,invli2n,logh,MPFR_RNDN);
  mpfr_set(den_bn,n_mpfr,MPFR_RNDN);
  mpfr_log(logn, n_mpfr, MPFR_RNDN);
  mpfr_mul(den_bn, den_bn, logn, MPFR_RNDN);
  mpfr_root(den_bn, den_bn, 4, MPFR_RNDN);
  mpfr_div(bn, num_bn, den_bn, MPFR_RNDN);
  mpfr_printf ("%.30RZf\n", bn);
}

int main(int argc, char* argv[]){
  mpz_t n;
  int Pplus, Log, Factors, Value, ShowBvalue;

  Pplus=Log=Factors=Value=ShowBvalue=0;

  int c;
  while((c=  getopt(argc, argv, "lpfvb")) != EOF) {
    switch (c) {
    case 'p':
      Pplus=1;
      break;
    case 'l':
      Log=1;
      break;
    case 'f':
      Factors=1;
      break;
    case 'v':
      Value=1;
      break;
    case 'b':
      ShowBvalue=1;
      break;
    }
  }

  mpz_init(n);
  mpz_set_str(n, argv[optind], 10);

  long p=atol(argv[optind+1]);
  long m=atol(argv[optind+2]);
  Hn hn(n,p,m);
  mpfr_printf("hn initialise  %.30RZ\n",hn.n_mpfr);
  
  //  printf("n= %ld p= %ld   m=%ld\n",n,p,m);
  //Gdelta res(p,m,maxdinit,verbose);

  if (Factors) {
    hn.gdelta.show_factors();
    }
  
  if (Log) {
        mpfr_printf("Compute Log h(n):\n");
	hn.compute_Log();
	mpfr_printf("Log h(n): %.26Re\n\n",hn.logh);
  }

  if (Pplus)
    cout << "Pplus = " << hn.gdelta.Pplus() << endl << endl;
  if (Factors) {

  }
  if (ShowBvalue) {
    hn.compute_bn();
  }
return 0;
}
