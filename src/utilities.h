#ifndef _utilities_h
#define _utilities_h
#include<iostream>
#include<vector>
#include<gmpxx.h>
#include<cmath>
using namespace std;
void erreur(string);

template< class T> void show(T x, size_t cnte) {
  cnte = min(cnte,x.size());
  cout << "{";
  for (unsigned int i=0; i< cnte-1; i++) {
    cout<< x[i] << ", ";
  }
  cout<< x[cnte-1] << "}\n";
}

template<class T> int find(T* table, T value, int imin, int imax);
#include"find.cpp"
#endif
