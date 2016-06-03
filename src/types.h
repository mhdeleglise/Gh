#ifndef _types_h
#define _types_h

#include<gmpxx.h>
#include<iostream>
#include<vector>
#include"utilities.h"

using namespace std;

struct kpkdelta{
  int k;
  int p;
  int delta;
  int sum;
  void display() {
    cout << "[" << this->k << ", " << this->p << ", " << this-> sum << ", " << this->delta << "]";
  }
};

kpkdelta pkd(long);

class GcombClass{
 public:
  GcombClass(int,long,int);
  void  display(int);
  int pk; // le pivot
  int Rhat;
  int k; // Le rang du pivot
  int P1; // Le premier des petits premiers
  int Pr; // Le plus grand des premiers
  int Increm;
  int Maxpgen;
  int Pmaxreel;
  int Maxm;
  vector<int> P;
  vector<mpq_class> Gvalues;
};

#endif
