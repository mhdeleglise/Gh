#include<gmpxx.h>
#include<cmath>
#include<vector>
#include"utilities.h"
#include"Gfunction.h"
#include<cmath>
#include<gmp.h>

void prevprime(mpz_t rop, mpz_t x);


int main(int argc, char* argv[]){
  if (argc != 4) {
    std::cout << "Usage: factor pk num den\n";
    exit(1);
  }
  //printf("%s  %s  %s\n", argv[1], argv[2], argv[3]);
  mpz_t pk, q, num, den;
  mpz_inits(pk,q,num, den, (mpz_ptr) 0);
  mpz_set_str(pk, argv[1],10);
  mpz_set_str(num, argv[2],10);
  mpz_set_str(den, argv[3],10);
  mpz_set(q, pk);
  
  int cnte=0;
  printf("Nk= 2 x 3 x 5 x .... x %s\n",argv[1]);
  

  gmp_printf("Num: ");
  while (mpz_cmp_si(num,1)) {
    while (!mpz_divisible_p(num, q)) {
      mpz_nextprime(q, q);
    }
    cnte+=1;
    mpz_divexact(num, num, q);
    //gmp_printf("Q= %.Zd  -> num = %.Zd\n",q,num); 
    gmp_printf("%.Zd ",q); 
  }
  gmp_printf("\n");
  //printf("Here we are cnte = %d\n",cnte);

    
  // Looking for the divisors Q1, Q2,     ,Qs  starting from p_[k+1], .....

  
  gmp_printf("Den: ");
  if (cnte==1) {
    printf("%s\n",argv[3]);
    return 0;
    }

  
  int i=0;
  while (i < cnte-1) {
    while (!mpz_divisible_p(den, pk)) {
      prevprime(pk, pk);
    }
    i+=1;
    mpz_divexact(den, den, pk);
    //gmp_printf("Q= %.Zd  -> den = %.Zd\n",pk,den);
    if (i < cnte-1)
      gmp_printf("%.Zd ",pk);
    else {
      gmp_printf("%.Zd ",pk);
      gmp_printf("%.Zd\n",den);
    }
  }

  return 0;
}
