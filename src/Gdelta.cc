#include<gmpxx.h>
#include<cmath>
#include<vector>
#include"numTheory.h"
#include"utilities.h"
#include"Gfunction.h"

//Calculateurs combinatoires globaux

Gdelta::Gdelta(long pk,long m,int maxdinit,int imp){
  init(pk,m,maxdinit,imp);
}

void Gdelta::init(long pk,long m,int maxdinit,int imp) {
  qlist.resize(0);
  dlist.resize(0);
  this->pk=pk;
  this->m = m;
  maxdinit=maxdinit-maxdinit%2;

  if (!is_prime(pk)) {
    erreur("Gpkmh, pk n'est pas premier\n");
  }
  if (imp >= 3)
    cout << "pk= " << pk << "    nprime= " << m << endl;

  long pkp1=next_prime(pk);
  mpq_class unsurpkp1(1,pkp1);

  if ((m>pkp1) || (m<0)) {
    erreur("Gpgmh, m is wrong");
  }

  if (m >= pkp1-2) {
    type= 0;
    ecart = m-pkp1+2;
    Gprov=pkp1;
    Gprov/=2;
    return;
  }
  
  if (m < pkp1-pk) {
    type= 1;
    ecart=m;
    Gprov=1;
    if (imp>=2) {
      cout << "type=1 m is smaller than pkp1-pk\n";
    }
    return;
  }

  long mm=m-m%2; // G(pk,m) = G(pk,mm)
  long pkp1mmm= pkp1-mm;
  if (is_prime(pkp1mmm)) {
    if (imp>=2)
      cout << "pkp1mmm= " << pkp1mmm << " is prime " << endl;
    type= 2;
    ecart= m%2;
    Gprov=pkp1;
    Gprov/=pkp1mmm;
    if (imp>=2)
      cout << "type=2 pkp1-mm is prime\n";
    return;
  }

  long pkp2=next_prime(pkp1);
  long maxdelta;
  if (mm%9 == 0)
    maxdelta= 2*(mm/9)-2;
  else
    maxdelta= 2*(mm/9);
  if (imp >= 6)
    cout << "Plus grand nombre pair < 2/mm/9 =  " << maxdelta << endl;
  // maxdelta is the greatest even number < 2m/9
  // 2*floor(1.275*lgpk*lgpk) is an empirical bound for delta_1(pk) defined in ZZZ
  double lgpk=std::log(pk);
  long delta1max=2*floor(1.275*lgpk*lgpk);
  maxdelta=min(maxdelta,delta1max);

  long maxd=min(pkp2-3,long(maxdinit));
  if (imp>=2)
    cout << "maxdinit= " <<  maxdinit <<  "  delta1max= " << delta1max << "  maxdelta= " << maxdelta << endl;
  long q=next_prime(pkp1mmm); // Rappel : ici pkp1mmm n'est pas premier
  if(imp>+2)
    cout << "q= " << q << endl;
  int oldd;
  long d=q-pkp1mmm; // smallest value for delta
  int delta_found=false;
  Gprov=1;
  while ((d <= maxdelta) && (not delta_found)) {
    qlist.push_back(q);
    dlist.push_back(d);
    if (d > maxd) {
      while(d>maxd)
	maxd=min(2*maxd,pkp2-3);
    }
    G1.init(pkp1,maxd,imp);
    mpq_class Gpr(G1[d]);
    Gpr*=pkp1;
    Gpr/=q;
    if (Gpr > Gprov) {
      Gprov=Gpr;
      dmax=d;
      qmax=q;
      cnteqmax=qlist.size();
    }
    if(imp>=2)
      cout << "d=" << d << "   G1[d] = " << G1[d] << "    1+d/pkp1 = "	\
	   <<  1+ d*unsurpkp1 << endl;;
    if (G1[d] >= 1+d*unsurpkp1) {
      delta_found = true;
      cntedmax=qlist.size();
      delta=d;
    }
    q = next_prime(q);
    oldd=d;
    d=q-pkp1mmm;
  }    
  // If we went out of the loop because delta is found, the Gdelta-method worked  
if (delta_found) {
    delta=oldd;
    // Sans le long double, une erreur dans le calcul de h(10^35), bien que le résultat soit vrai par chance
    // la valeur de double de pkp1 perd les 3 derniers digits
    cout.precision(25);
    //cout << "double pkp1 = " << (double)pkp1 << endl;
    long qhat=(long double)(pkp1)*pkp2*(pkp1mmm+delta)/(pkp1+delta)/(pkp1-3*delta/2);
    if(imp>=4) {
      cout << "imp="<< imp << endl;
      cout << "pkp1  pkp2 " << pkp1 << " " << pkp2 << endl;
      cout << "pkp1mmm+delta=" << pkp1mmm+delta << endl;
      cout << "pkp1+delta=" << pkp1+delta << endl;
      cout << "pkp1-3*delta2 =" << pkp1-3*delta/2 << endl;
      cout << "delta is found, quat= " << qhat << endl;
    }
    long maxiq=previous_prime(qhat+1);
    maxid=maxiq-pkp1mmm;
    if (maxid > maxd) {
      cout << "maxd incremente de " << maxd << " a " << maxid << endl;
      maxd=maxid;
      G1.init(pkp1,maxd,0);
    }
    if (imp>=2)
      cout << "Calcul de G(pk,mm) = max(pkp1/q*G(pkp1,q-pk1mm), dernieres valeurs de q \n";
    while (q <= qhat) {
      qlist.push_back(q);
      dlist.push_back(d);
      mpq_class Gpr=pkp1;
      Gpr/=q;
      Gpr*=G1[d];
      if (Gpr > Gprov)
	{
	  Gprov=Gpr;
	  dmax=d;
	  qmax=q;
	  cnteqmax=qlist.size();
	  if (imp>=2)
	    cout << "Gprov > Gpr nouvelle valeur de dmax = d " << dmax << endl;
	}
      q=next_prime(q);
      d=q-pkp1mmm;
    }
    // "Calcul de ecart\n";
    d=dmax-2;
    while ((d >= 0) && (G1[d] == G1[dmax])){
      d -= 2;
    }
    ecart=dmax-d-2+m%2;
    type = 3;
  }
// If we didn't find delta
 else {
   G2.init(pk,mm,imp);
   Gprov=G2[m];
   d=mm-2;
   while (G2[d] == Gprov)
     d -= 2;
   ecart = mm-d-2+m%2;
   type=4;
   if (imp > 1) {
     if (mm > 9*delta1max/2)
       erreur("The Gdelta-method has failed");
     else
       cout << "The Gdelta-method has failed but m is small" << endl;
   }
 }
 return; 
}


void Gdelta::display(int verbose) {
  if (!verbose) {
  cout << "[" << Gprov << ", " << ecart << ", " << type;
  switch(type) {
  case 3:
    cout << ", " << delta << ", " << maxid  << ", " << qlist.size();
    break;
  case 4:
    cout << ", " << m;
  }
  cout << "]" << endl;
  cout.precision(10);
  }
  else {
    cout << "Gprov= " << Gprov << endl;
    cout << "Ecart " << ecart;
    cout << "  type " << type;
    cout << "  delta " << delta;
    cout << "  maxid= " << maxid;
    cout << "  cnteq = " << qlist.size() << endl;;
    //cout << "  cntedmax= " << cntedmax << endl;
    cout << "dlist="; show(dlist,cntedmax);
    cout << "qlist=";show(qlist,qlist.size());
  }
}

void Gdelta::show_pk_Gp() {
  cout << pk << " " << Gprov.get_num() << " " << Gprov.get_den() << endl;
}

long double Gdelta::log() {
  mpz_class num = Gprov.get_num();
  mpz_class den = Gprov.get_den();
  return std::log((long double)num.get_d())-std::log((long double)den.get_d());
}

void sliceGdelta(long pk, long sk) {
  long pkp1=next_prime(pk);
  int m,lastm=pkp1-2;
  Gdelta res;
  for(m=0; m < lastm; m+= 2) {
    res.init(pk,m,500,1);
    cout << "m= " << m << "  x= " << sk+m << " Gdelta= ";
    res.display();
    cout << endl;
  }
  m=lastm;
  res.init(pk,m,500,1);
    cout << "m= " << m << "  x= " << sk+m << " Gdelta= ";
    res.display();
    cout << endl;
}

void Gdelta::show_log() {
  cout << "Log   = " << log() << endl << endl;
}

void Gdelta::show_frac() {
  cout << "Gfrac = " << Gprov << endl << endl;;
}

long Gdelta::Pplus() {
  mpz_class numer=Gprov.get_num();
  long pmax=pk;
  mpz_t q;
  if (numer > 1) {
    mpz_init_set_si(q,pk);
    mpz_nextprime(q,q);
    while (numer > 1)
      {
	while (! mpz_divisible_p(numer.get_mpz_t(), q)) {
	  mpz_nextprime(q,q);
	}
	mpz_divexact(numer.get_mpz_t(), numer.get_mpz_t(), q);
      }
    pmax=mpz_get_si(q);
  }
  return pmax;
}

void prevprime(mpz_t rop, mpz_t x);

void Gdelta::show_factors() {
  int cnte=0;
  mpz_class numer=Gprov.get_num();
  mpz_t q;
  if (numer > 1) {
    cout << "Facteurs premiers du numerateur   : ";
    mpz_init_set_si(q,pk);
    mpz_nextprime(q,q);
    while (numer > 1)
      {
	cnte += 1;
	while (! mpz_divisible_p(numer.get_mpz_t(), q)) {
	  mpz_nextprime(q,q);
	}
	gmp_printf("%.Zd ", q);
	mpz_divexact(numer.get_mpz_t(), numer.get_mpz_t(), q);
      }

    cout << endl;
    cout << "Facteurs premiers du denominateur : ";
    int i=0;
    mpz_class denom=Gprov.get_den();
    mpz_init_set_si(q,pk);

  if (cnte==1) {
    gmp_printf("%.Zd\n\n",denom.get_mpz_t());
    return;
    }

    
    while (i < cnte-1) {
      while (!mpz_divisible_p(denom.get_mpz_t(), q)) {
	prevprime(q,q);
      }
      i+=1;
      mpz_divexact(denom.get_mpz_t(), denom.get_mpz_t(), q);
      if (i < cnte-1)
	gmp_printf("%.Zd ",q);
      else {
	gmp_printf("%.Zd ",q);
	gmp_printf("%.Zd\n",denom.get_mpz_t());
      }
    }
  }
  cout << endl;
}
