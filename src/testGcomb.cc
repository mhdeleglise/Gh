#include<iostream>
#include<iomanip>
#include"Gfunction.h"
//#include"types.h"
using namespace std;
long next_prime(long);
int is_prime(long);


int main(int argc, char* argv[]) {
  if ((argc==1) || (argc > 4))
    {
      cout << "Usage testGcomb pk mmax verbose\n";
      cout << "Calcule G(pk,m), 0 <= m <= mmax, par la méthode combinatoire\n";
      cout << "pk est premier\n";
      cout << "verbose, facultatif, est un entier de [0..8] (0 par défaut)\n";
      exit(0);
    }
    
  long p = atol(argv[1]);
  int maxm = atoi(argv[2]);
  int imp;
  if (argc<4)
    imp=0;
  else
    imp = atoi(argv[3]);
  cout << "Appel de Gcomb avec p= " <<p<< "  maxm=" << maxm << "  imp=" << imp  << endl;
  GcombData data(p,maxm,imp);
  data.display(imp);
}
