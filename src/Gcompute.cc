#include<gmpxx.h>
#include<cmath>
#include<vector>
#include<unistd.h>
#include"utilities.h"
#include"Gfunction.h"
#include"Nk.h"
int main(int argc, char* argv[]){

  int Pplus, Log, Factors, Value;
  const int maxdinit=500;
  const int verbose=0;

  Pplus=Log=Factors=Value=0;

  int c;
  while((c=  getopt(argc, argv, "lpfv")) != EOF) {
    switch (c) {
    case 'p':
      Pplus=1;
      break;
    case 'l':
      Log=1;
      break;
    case 'f':
      Factors=1;
      break;
    case 'v':
      Value=1;
      break;
    }
  }
  
  long p=atol(argv[optind]);
  long m=atol(argv[optind+1]);Gdelta res(p,m,maxdinit,verbose);

  if (Factors)
    res.show_factors();

  if (Log)
    res.showLog();
  if (Pplus)
    cout << "Pplus = " << res.Pplus() << endl << endl;
  if (Factors) {

  }
  if (Value) {
    res.show_value();
  }
return 0;
}
