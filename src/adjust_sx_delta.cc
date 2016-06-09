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
  cout << "Adjust x= " << x << endl;
  mpz_t deltax;
  mpz_init(deltax);
  mpz_set_str(deltax, argv[index+1],10);
  gmp_printf("x= %ld delta= %Zd\n",x,deltax);
  long64 maxprime = min(sqrt((double)(2*x))+10000, 2000000000.0);
  long64 p=0;
  presieved_primes::init_prime_table(maxprime,2);
  prime_generator pg(10000000, x);

  mpz_t sum, delta;
  mpz_init(sum);
  mpz_init(delta);
  long64  cnte=1;
  long pj1=pg.next_prime(x);
  cout <<  "pj1= " << pj1 << endl;
  gmp_printf("x= %ld   deltax= %.Zd\n",x,deltax);
  if (mpz_cmp_si(deltax, 0) >= 0) {
    long pj;
    while (mpz_cmp_si(deltax, pj1) >= 0) {
      mpz_sub_ui(deltax, deltax, pj1);
      cnte+=1;
      pj=pj1;
      pj1=pg.next_prime();
    }
    //gmp_printf("cnte= %d pk= %ld deltax= %.Zd\n",cnte,pj,deltax);
    if (verbose)
      gmp_printf("cnte= %d pk= %ld deltax= %.Zd\n",cnte,pj,deltax);
    else
      gmp_printf("%ld %.Zd\n",pj,deltax);
    return 0;
    }
  else {
    p= pg.prev_prime(x);
    mpz_add_ui(deltax, deltax, p);
    while (mpz_cmp_si(deltax, 0) < 0) {
      p=pg.prev_prime();
      cnte++;
      mpz_add_ui(deltax, deltax, p);
    }
    mpz_set(delta, deltax);
    p=pg.prev_prime();
    if (verbose)
      gmp_printf("cnte= %d pk= %ld delta= %.Zd\n",cnte,p,delta);
    else
      gmp_printf("%ld %.Zd\n",p,delta);
    return 0;
  }
}
