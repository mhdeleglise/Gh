#include <iostream>
#include <sstream>
#include <iomanip>
#include<gmp.h>


using namespace std;


int main(int argc, char* argv[]) {
  mpz_t x1;
  mpz_t x2;
  mpz_init(x1);
  mpz_init(x2);

  mpz_set_str(x1, argv[1], 10);
  mpz_set_str(x2, argv[2], 10);
  mpz_sub(x1,x1,x2);

  gmp_printf ("%Zd\n", x1);
  return 0;
}
