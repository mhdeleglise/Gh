#include<iostream>
#include<iomanip>
#include"Gfunction.h"
#include"types.h"
using namespace std;
long next_prime(long);
int is_prime(long);

GcombData::GcombData(long pivot, int maxm, int imp=0) {
  init(pivot,maxm,imp);
}

GcombData::GcombData() {
  P.resize(0);
  Sp.resize(0);
  Gvalues.resize(0);
  Maxm=0;
  K=0;
  Pk=0;
  Pkp1=0;
}

void GcombData::init(long pivot, int maxm, int imp=0) {
  if ((Pk==pivot) && (maxm <= Maxm)) {
    //cout << "Gcombdata::init Rien à faire\n";
    return;
}
  this->rhat=Gfun::rhat_default;
  computation_done=0;
  P.clear();
  Sp.clear();
  Gvalues.clear();
  HH.resize(Gfun::rhat_default+1);

  if (is_prime(pivot) == 0) {
    cout << "pivot = " << pivot << endl;
    erreur("Error in Gcombdata constructor pivot is not prime\n");
  }
  this->Pk = pivot;
  this->Maxm=maxm;
  this->Pkp1 = next_prime(Pk);


  if (Maxm > Pkp1-3) {
    cout << "Maxm = " << Maxm  << " bigger than  " << Pkp1-3 << endl;
    erreur("is too big ");
  }
  if (Maxm < Pkp1-Pk) {
    int m2;
    Gvalues.resize(Maxm/2+1);
    for (m2=0; m2 <=Maxm/2; m2++)
      Gvalues[m2]=1;
    K=0;
    return;
  }

  make_primes();
  do {
    computes(imp);
    }
  while (!computation_done);
}

void GcombData::make_primes() {
  this-> P1 = next_prime(Pkp1-Maxm-1);
  long Pr = P1;

  Sp.push_back(P1);
  P.push_back(P1);
  while (Pr <= Pk) {
    P.push_back(Pr);
    Sp.push_back(Pr+Sp[Sp.size()-1]);
    Pr = next_prime(Pr);
  }
  K = P.size()-1;
  //cout << "Les grands nombres premiers\n";
  for (int r=0; r <= rhat; r++) {
    P.push_back(Pr);
    Sp.push_back(Pr+Sp[K+r]);
    Pr = next_prime(Pr);
  }
} 

void GcombData::one_more_big_prime(int verbose) {
  if (verbose>=4)
    cout << "In one more big prime \n";
  long Pr= P[P.size()-1];
  if (verbose>=4)
    cout << "P.size=" << P.size() << " Dernier elet de P=" << Pr << endl;
  Pr=next_prime(Pr);
  P.push_back(Pr);
  Sp.push_back(Pr+Sp[K+rhat+1]);
  rhat=rhat+1;
}


void GcombData::computes(int imp){
  if (imp >=2)
    cout << "In computes avec rhat = " << rhat << endl;
  mpz_class pro = 1;
  int r,j;
  long m2;
  for (r=K+1; r <= K+rhat ; r++) {
    pro = pro*P[r];
  }
  if (rhat > Gfun::rhat_default) {
    HH.resize(rhat+1);
}
  for (j=0; j <= this->rhat; j++) {
    HH[j].assign(Maxm/2+1,0);
  }

  for (m2=0; m2 <= Maxm/2; m2++) {
    HH[0][m2] = 1;
  }
  
  int j1;
  long pkjr;
  //   LA BOUCLE PRINCIPALE   
  if (imp >=5)
    cout << "Dans la boucle principale pour r de 1 a K+rhat = " << K+rhat << endl;
  for (r=1; r <= K+rhat ; r++) {
    if (imp >= 6)
      cout <<"r = " << r << endl;
    for (j=min(K+rhat-K,r); j >= max(1,r-K);  j--) {
      j1 = j-1; 
      long mj1 = mj(j,K,r-1);
      pkjr = (P[K+j]-P[r])/2;
      if (imp >= 8)
	cout << "       j=" << j << "  mj1 = " << mj1 << " pkjr = " << pkjr << endl;
      long upbound;
      if (mj1==Gfun::infinity) {
	upbound = Maxm/2;
      }
      else
	upbound = min(Maxm/2,mj1-1);
      if (imp >= 6)
	cout << "          Boucle 1" << endl;
      for (m2=mj(j,K,r); m2 <= upbound; m2++) {
	//	cout << "m2= " << m2 << " j1= " << j1 << "  m2-pkjr=" << m2-pkjr << endl;
	HH[j][m2] = P[r]*HH[j1][m2-pkjr];
      }
      if (imp >= 6)
	cout << "          Boucle 2" << endl;
      for (m2=mj1; m2 <= Maxm/2; m2++) {
	if (HH[j][m2] > P[r]*HH[j1][m2-pkjr]) {
	  HH[j][m2] = P[r]*HH[j1][m2-pkjr];
	}
      }
      
      if (imp >= 6)
	cout << endl;
    }
  }
  // RECOLLECTION DES VALEURS DE GG et determination of maxPgen
  if (imp >= 5) {
    cout << "\nRecollection des valeurs de G et calculde maxpgen" << endl;
  }

  Maxpgen = 0;
  long maxP = 0;
  Gvalues.resize(Maxm/2+1);
  for (m2 = 0; m2 <= Maxm/2; m2++) {
    if (imp >= 6)
      cout<< "m= " << setw(4) << 2*m2 << "   ";
    Gvalues[m2] = mpq_class(pro,HH[this->rhat][m2]);
    Gvalues[m2].canonicalize();
    if (imp >= 8)
      cout <<  "    HH[rhat][m2] = " << HH[this->rhat][m2] << endl;
    if (imp >= 6)
      cout << " Gvalues[m] = " << Gvalues[m2] << endl;
    if (Gvalues[m2] == 1) 
      maxP = 2*m2+P[K];
    else{
      mpq_class tmp = Gvalues[m2]/(Gvalues[m2]-1);
      if (2*m2+P[K] <= Gvalues[m2]*2*m2/(Gvalues[m2]-1)) { 
	maxP = 2*m2+P[K];
      }
      else
	{
	  mpq_class tmp = 2*m2*Gvalues[m2]/(Gvalues[m2]-1);
	  mpz_class num = tmp.get_num();
	  mpz_class den = tmp.get_den();
	  mpz_class quot(1);
	  mpz_fdiv_q(quot.get_mpz_t(),num.get_mpz_t(),den.get_mpz_t());
	  maxP = quot.get_si();
	}
    }
    if (maxP > Maxpgen) {
      // cout<<  "Maxpgen passe de " << maxPgen << " a    " << maxP << endl;
      Maxpgen = maxP;
    }
  }
  if (imp >= 5) {
    cout << "maxpgen = " << Maxpgen << "     P[K+rhat+1] = " << P[K+rhat+1] << endl;
  }

  if (Maxpgen >= P[K+rhat+1]) {
    if (imp >= 2)
      cout << "rhat= " << rhat << " est trop petit\n";  
    if (imp >= 4)
      cout << "maxPgen = " << Maxpgen  << "   >= P[K+rhat+1] = " << P[K+rhat+1] << endl; 
    one_more_big_prime(imp);
    return; // this call of computes id finished. 
    //  Computes() is going to start again with one more prime.
  }    
  else {
    computation_done=1;
    int ilast = find<long>(&P[K],Maxpgen,1,this->rhat+1)-1;
    this->RealIncrem=ilast;
  
    // We measure the largest prime which has been used
    int rmes=K+1;
    int mmes=0;
    for (m2=0; m2 <= Maxm/2; m2++) {
      mpz_class num=(Gvalues[m2]).get_num();
      for (r=rmes+1; r < K+rhat; r++)
	if (num % P[r]==0) {
	  rmes=r; 
	  mmes=2*m2;
	}
    }
    this->Increm = rmes-K;
    this->Pmax = P[rmes];
    this->mPmax=mmes;
    Sp.clear();
    P.clear();
  }
}


void GcombData::display(int verbose) {
  if (verbose >= 3) {
    cout << "\nG(pk,m) values for pk= " << Pk << "  0 <= m <= " << Maxm << ":\n";
    int m2;
    for (m2 = 0; m2 <= Maxm/2; m2++) {
      cout << 2*m2 << ":" << Gvalues[m2] << " ";
      cout << endl;
    }
    cout << endl;
  }
  
  int j;
  if (K > 0) {
    cout << "Pk= " << Pk << "  Maxm= " << Maxm << "  rhat=" << rhat << endl;
  }
  if (K > 1) {
  cout << "K= " << K << "  K+rhat= " << K+rhat << endl;
  cout << "P1= " << P[1] << "   P[K+rhat] = " << P[K+rhat] << "   P[K+rhat+1] = " << P[K+rhat+1] << endl;
  cout << "Les grands premiers" << endl;
  cout << "{";
  for (j=K+1; j < K+rhat; j++) 
    cout << P[j]  <<", ";
  cout << P[K+rhat] << "}" << endl;
  cout << "Increm     = " << this->Increm << endl;
  cout << "RealIncrem = " << this->RealIncrem << endl;
  cout << "Pmax       = " << this->Pmax << endl;
  cout << "mPmax      = " << this->mPmax << endl;
  cout << "maxPgen    = " << this->Maxpgen << endl;
  }
  else
    cout << "K=0, que des 1\n";
}

