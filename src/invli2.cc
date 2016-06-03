#include<iostream>
#include<gmp.h>
#include<mpfr.h>
#include<math.h>
#include <stdlib.h>     /* atof */
using namespace std;

// $1=n
// Computes x=sqrt(Li^(-1)(n))
// ie Li(x^2)=n

void invli(mpfr_t res, mpfr_t x, double precision);
  
int main(int argc, char* argv[]) {
  mpfr_t n,x2,x1;
  mpfr_inits2(256,n,x1,x2, (mpfr_ptr) 0);
  mpfr_set_str (n, argv[1], 10, MPFR_RNDN);

  double precision= (argc > 2) ? atof(argv[2]) : 1e-12;
  invli(x2,n,precision);
  mpfr_sqrt(x1,x2,MPFR_RNDN);
    mpfr_printf ("%.30RZf\n", x1);
}
