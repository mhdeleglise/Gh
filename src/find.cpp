#include<iostream>
#include"utilities.h"

template<class T> int find(T* table, T value, int imin, int imax) {
  if ((imax < imin) || (table[imin] > value)|| (imin < 0)) {
    erreur("Donnes incorrectes pour find");
    cout << "imin " << imin << "  table[imin] = " << table[imin] << "  imax " << imax << "  table[imax] " <<  table[imax] << endl;
    cout << "Value = " << value << endl;
  }
  if (value >= table[imax])
    return imax;

  int a = imin;
  int b = imax;
  while (b-a > 1) {
    int c = (a+b)/2;
    if (table[c] <= value)
      a = c;
    else
      b = c;
  }
  return a;
} 
