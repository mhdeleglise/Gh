#include<iostream>
#include<iomanip>
#include<math.h>
#include<mpfr.h>

char const* gamma_euler="0.5772156649015328606065120900824024310421593359399235988057672348848677267777";

using namespace std;

inline double max(double a, double b) {
  return (a > b) ? a : b;
}


void li(mpfr_t res, mpfr_t x, double precision) {
  mpfr_t local_res, a, b, ab, som1, som2, un_surn, delta, gamma, lgx, loglogx, sqrtx;
  mpfr_inits2(128, local_res, a, b, ab, som1, som2, un_surn, delta, gamma, lgx, loglogx, sqrtx, (mpfr_ptr) 0);

  mpfr_log(lgx, x, MPFR_RNDN);
  mpfr_log(loglogx, lgx, MPFR_RNDN); 
  mpfr_sqrt(sqrtx, x, MPFR_RNDN);
  mpfr_set_str(gamma, gamma_euler, 10, MPFR_RNDN);
  mpfr_add(local_res, gamma, loglogx, MPFR_RNDN);

  /*
  mpfr_printf("lgx       = %.30Rf\n",lgx);
  mpfr_printf("loglogx   = %.30Rf\n",loglogx);
  mpfr_printf("sqrtx     = %.30Rf\n",sqrtx);
  mpfr_printf("gamma     = %.30Rf\n",gamma);
  mpfr_printf("local_res initialise : local_res = gamma + loglogx = %.30Rf\n",local_res);
  */

  double delta_bound=precision/mpfr_get_d(sqrtx,MPFR_RNDN);

  static int n1= max(8,ceil(5*mpfr_get_d(lgx,MPFR_RNDN)/9));
  
  if (n1 % 2 == 0)
    n1 += 1;
  //cout << "n1 = " << n1 << endl;

  mpfr_set(a, lgx, MPFR_RNDN);
  mpfr_set_si(b, 1, MPFR_RNDN);
  mpfr_set(som1, lgx, MPFR_RNDN);
  mpfr_mul(a, a, lgx, MPFR_RNDN);
  mpfr_div_si(a, a, 4, MPFR_RNDN);

  mpfr_set(som2,som1, MPFR_RNDN);
  mpfr_mul(ab, a, b, MPFR_RNDN);
  mpfr_sub(som2, som2, ab, MPFR_RNDN);


  int n=3;
  for (; n <= n1; n+= 2) {
    //a *= lgx/(2*n);
    mpfr_mul(a, a, lgx, MPFR_RNDN);
    mpfr_div_si(a, a, 2*n, MPFR_RNDN);

    //b += 1.0/n;
    mpfr_set_si(un_surn, 1, MPFR_RNDN);
    mpfr_div_si(un_surn, un_surn, n, MPFR_RNDN);
    mpfr_add(b, b, un_surn, MPFR_RNDN);

    //som1 = som2+a*b;
    mpfr_mul(ab, a, b, MPFR_RNDN);
    mpfr_add(som1, som2, ab, MPFR_RNDN);

    //a *= lgx/(2*n+2);
    mpfr_mul(a, a, lgx, MPFR_RNDN);
    mpfr_div_si(a, a, 2*n+2, MPFR_RNDN);

    
    //som2 = som1-a*b;
    mpfr_mul(ab, a, b, MPFR_RNDN);
    mpfr_sub(som2, som1, ab, MPFR_RNDN);  
  }

  n -= 2;
  //printf("HERE\n");
  mpfr_sub(delta, som1, som2, MPFR_RNDN);
  //mpfr_printf("delta= %.30Rf\n", delta);
  
  while (mpfr_get_d(delta, MPFR_RNDN) > delta_bound) {
    n+= 2;
    //a *= lgx/(2*n);
    mpfr_mul(a, a, lgx, MPFR_RNDN);
    mpfr_div_si(a, a, 2*n, MPFR_RNDN);

    //b += 1.0/n;
    mpfr_set_si(un_surn, 1, MPFR_RNDN);
    mpfr_div_si(un_surn, un_surn, n, MPFR_RNDN);
    mpfr_add(b, b, un_surn, MPFR_RNDN);

    //som1 = som2+a*b;
    mpfr_mul(ab, a, b, MPFR_RNDN);
    mpfr_add(som1, som2, ab, MPFR_RNDN);

    //a *= lgx/(2*n+2);
    mpfr_mul(a, a, lgx, MPFR_RNDN);
    mpfr_div_si(a, a, 2*n+2, MPFR_RNDN);

    
    //som2 = som1-a*b;
    mpfr_mul(ab, a, b, MPFR_RNDN);
    mpfr_sub(som2, som1, ab, MPFR_RNDN);  
    mpfr_sub(delta, som1, som2, MPFR_RNDN);
  }
  mpfr_mul(ab, sqrtx, som1, MPFR_RNDN);
  mpfr_add(local_res, local_res, ab, MPFR_RNDN);
  mpfr_set(res, local_res, MPFR_RNDN);
}


void invli(mpfr_t res, mpfr_t x, double prec) {
  mpfr_t y, logy, liy, lgx, tn, logtn, delta;
  mpfr_inits2(128, y, logy, liy, lgx, tn, logtn, delta, (mpfr_ptr) 0);
  mpfr_log(lgx, x, MPFR_RNDN);
  mpfr_mul(y, x, lgx, MPFR_RNDN);
  li(liy, y, prec);
  
  mpfr_log(logy, y, MPFR_RNDN);
  mpfr_sub(delta, x, liy, MPFR_RNDN); 
  while (mpfr_get_d(delta, MPFR_RNDN) > prec) {
    mpfr_mul(delta, delta, logy, MPFR_RNDN);
    mpfr_add(y, y, delta, MPFR_RNDN);
    li(liy, y, prec);
    mpfr_log(logy, y, MPFR_RNDN);
    mpfr_sub(delta, x, liy, MPFR_RNDN); 
  }
  mpfr_set(res, y, MPFR_RNDN);
}
