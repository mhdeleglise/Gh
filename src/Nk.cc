#include<stdio.h>
#include<iostream>
#include<gmp.h>
#include<mylib.h>



void Nk_compute(mpz_t Nk, long pk) {
  mpz_t res;
  mpz_init_set_si(res,6);
  presieved_primes::init_prime_table(1000000,2);
  prime_generator pg(10000000);
  long p=pg.next_prime();
  while (p <=pk) {
    mpz_mul_si(res,res,p);
    p=pg.next_prime();
  }
  mpz_set(Nk,res);
}


