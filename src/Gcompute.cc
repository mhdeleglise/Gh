#include<gmpxx.h>
#include<cmath>
#include<vector>
#include<unistd.h>
#include"utilities.h"
#include"Gfunction.h"

int main(int argc, char* argv[]){

  int Verbose, Pplus, Log, Factors;
  const int maxdinit=500;
  const int verbose=0;

  Pplus=Log=Factors=0;

  int c;
  while((c=  getopt(argc, argv, "vlpf")) != EOF) {
    switch (c) {
    case 'v':
      Verbose=1;
      break;
    case 'p':
      Pplus=1;
      break;
    case 'l':
      Log=1;
      break;
    case 'f':
      Factors=1;
    }
  }
  
  long p=atol(argv[optind]);
  long m=atol(argv[optind+1]);
  Gdelta res(p,m,maxdinit,verbose);
  cout << endl;
  res.show_frac();
  if (Log)
    res.show_log();
  if (Pplus)
    cout << "Pplus = " << res.Pplus() << endl << endl;
  if (Factors) {
    res.show_factors();
  }
  return 0;
}
