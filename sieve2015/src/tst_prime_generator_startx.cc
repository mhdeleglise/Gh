#include<mylib.h>
#include<cmath>
#include"prime_generator.h"
#include<iostream>


int main(int argc, char* argv[]) {
  long x = (argc > 1) ? atol(argv[1]) : 100;
  long wsize = (argc > 2) ? atol(argv[2]) : 1000;

  presieved_primes::init_prime_table(1000000,2);
  presieved_primes::display();

  cout << "Creation d'un prime generator de taille " << wsize << "   startx = " << x << endl;
  prime_generator pg(wsize, x);

  long p=pg.next_prime();
  cout << "p= " << p << endl;

  if (p != 101)
    exit(1);

  return 0;
}
