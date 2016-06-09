#include<iostream>
#include<iomanip>
#include<gmp.h>
using namespace std;

void prev_prime(mpz_t x, char const *s ) {
  mpz_set_str(x, s, 10);

  if (mpz_cmp_si(x,2) == 0) {
    return;
  }
  


  if (mpz_divisible_2exp_p (x, 1))
    mpz_sub_ui(x, x, 1);

  while (! mpz_probab_prime_p (x, 100)) {
    mpz_sub_ui(x, x, 2);
  }
}

int main(int argc, char* argv[1]) {
  mpz_t x;
  mpz_init(x);
  prev_prime(x, argv[1]);
  gmp_printf("%Zd\n",x);
}
