#include<gmpxx.h>
#include<cmath>
#include<vector>
#include"utilities.h"
#include"Gfunction.h"

int main(int argc, char* argv[]){
  int verbose=0;
  int maxdinit=500;
  if ((argc==1) || (argc > 5)) {
    cout << "Usage testGdelta pk m maxdinit verbose\n";
    cout << "Computes Gdelta(pk,m)\n";
    cout << "default value for maxdinit=500, default value for verbose=0\n";
    exit(0);
  }
  long p=atol(argv[1]);
  long m=atol(argv[2]);
  if (argc>=4)
    verbose=atoi(argv[3]);
  if (argc>=5)
    maxdinit=atoi(argv[4]);

  Gdelta res(p,m,maxdinit,verbose);
  //res.show_pk_Gp();
  return 0;
}
