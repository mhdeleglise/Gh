#include<gmp.h>

void prevprime(mpz_t rop, mpz_t x) {
  if (mpz_divisible_2exp_p (x, 1))
    mpz_sub_ui(x, x, 1);
  else
    mpz_sub_ui(x, x, 2);
  mpz_set(rop,x);
  while (! mpz_probab_prime_p (rop, 100)) {
    mpz_sub_ui(rop, rop, 2);
  }
}
