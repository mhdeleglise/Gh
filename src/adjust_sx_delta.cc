// -*- C++ -*-
#include<mylib.h>
#include<math.h>
#include<stdlib.h>
#include<unistd.h>
#include<gmp.h>

int main(int argc, char * argv[]) {
  int verbose=0;
  int index;
  
  if (argc < 3) {
    cout << "Usage adjust_sx_delta [-v] x delta\n";
    exit(1);
  }
  if ((argv[1][0]=='-') && (argv[1][1]=='v')) {
    verbose=1;
    index=2;
  }
  else
    index=1;
  
  long64 x=atol(argv[index]);
  mpz_t sum_x;
  mpz_init(sum_x);
  mpz_set_str(sum_x, argv[index+1],10);
  long64 maxprime = min(sqrt((double)(2*x))+10000, 2000000000.0);
  long64 p=0;
  long64 lastp=p;
  presieved_primes::init_prime_table(maxprime,2);
  prime_generator pg(100000000, x);

  mpz_t sum, delta;
  mpz_init(sum);
  mpz_init(delta);
  long64  cnte=1;

  if (mpz_cmp_si(sum_x, 0) >= 0) {
    while (mpz_cmp(sum, sum_x) < 0) {
      lastp=p;
      p=pg.next_prime();
      mpz_add_ui(sum,sum,p);
      cnte++;
    }
    cnte--;
    mpz_sub_ui(sum,sum,p);
    mpz_sub(delta, sum_x, sum);  
    p=lastp;
    if (verbose)
      gmp_printf("cnte= %d p= %ld delta= %.Zd\n",cnte,p,delta);
    else
      gmp_printf("%ld %.Zd\n",p,delta);
    return 0;
  }
  else {
    p= pg.prev_prime(x);
    mpz_add_ui(sum_x, sum_x, p);
    while (mpz_cmp_si(sum_x, 0) < 0) {
      p=pg.prev_prime();
      cnte++;
      mpz_add_ui(sum_x, sum_x, p);
    }
    mpz_set(delta, sum_x);
    p=pg.prev_prime();
    if (verbose)
      gmp_printf("cnte= %d pk= %ld delta= %.Zd\n",cnte,p,delta);
    else
      gmp_printf("%ld %.Zd\n",p,delta);
    return 0;
  }
}
