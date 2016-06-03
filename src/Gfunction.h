#ifndef _Gfunction_h
#define _Gfunction_h
#include<iostream>
#include<vector>
#include<gmpxx.h>
#include<limits.h>
#include"utilities.h"
long next_prime(long);
long previous_prime(long);
int is_prime(long);

using namespace std;

namespace Gfun{
const int rhat_default = 10; 
const long infinity = LONG_MAX;
}

long next_prime(long);
int is_prime(long);

class GcombData{
 public:
  int rhat;
  long Pk;
  long P1;
  long Pkp1;
  long Maxm;
  int K;
  int Rhat;
  long Maxpgen; 
  long Pmax; // Le plus grand facteur premier de G(pk,mp pour 0 <= m <= Mmax
  int mPmax; // Le premier m pour laquelle G(pk,m) est divisible par Pmax
  int Increm;
  int RealIncrem;
  vector<mpq_class> Gvalues;
  void init(long pk, int maxm, int imp);
  mpq_class& operator[] (unsigned i) {return Gvalues[i/2];}
  GcombData(long pivot, int maxm, int imp);
  GcombData();
  void computes(int verbose=0);
  void display(int verbose=0);
  vector<vector<mpz_class> > HH;
  int computation_done;

 private:
  vector<long> P;
  //  vector<long> Sp;
  vector<mpz_class> Sp;
  // This procedure calculates m_j(P_r) defined in the article DNZBordeaux Equation (103)
  inline long mj(int j,int k,int r);
  void make_primes();
  void one_more_big_prime(int verbose);
};

class Gdelta{
  long pk;
  long m;
  int type;
  long dmax;
  long qmax;
  int maxid;
  long cnteqmax;
  long cntedmax;
  int delta;
  mpq_class Gprov;  // !!!!! Si Gprov est public cela plante
  vector<long> qlist;
  vector<long> dlist;
  GcombData G1;
  GcombData G2;
 public:
  long ecart;
  mpq_class value() {return Gprov;}
  long double log();
  Gdelta() {};
  Gdelta(long pk,long m,int maxdinit,int imp);
  void init(long pk,long m,int maxdinit,int imp);
  void display(int verbose=0);
  void show_pk_Gp();
};

inline long GcombData::mj(int j,int k,int r) {
  //  long a;
  mpz_class a;
  if (j>r)
    return Gfun::infinity;

  a = (Sp[k+j]-Sp[k]-Sp[r]+Sp[r-j])/2;
  if (a < 0) {
    cout << "in mj avec j,k,r= " << j << " " << k << " "  << endl;
    cout << "a = " << a << endl;
    cout << "k+j = " << k+j << "  Sp[k+j] = " << Sp[k+j] << endl;
    cout << "k  = " << k    << "  Sp[k  ] = " << Sp[k] << endl;
    cout << "r  = " << r    << "  Sp[r  ] = " << Sp[r] << endl;
    cout << "r-j = " << r-j << "  Sp[r-j] = " << Sp[r-j] << endl;
    erreur("Error mj a < 0");
    return 0;
  }
  else
    if (mpz_fits_slong_p(a.get_mpz_t()))
      return mpz_get_si(a.get_mpz_t());
    else {
      erreur("Dans Gcombdata::mj a=mj(j,k,r) est trop grand pour être convert en longint");
      return 0;
    }
}

#endif
