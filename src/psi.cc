#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <mpfr.h>
using namespace std;

mpfr_t _tmp1;
mpfr_t _big_sumlog;
mpfr_t _small_sumlog;
mpfr_t _sumlog;
mpfr_t _mcarre;
mpfr_t _ncarre;
mpfr_t _n;
mpfr_t _m;
mpfr_t _n1;
mpfr_t _m1;
mpfr_t _logn;
mpfr_t _logm;
mpfr_t _den;
mpfr_t _ti;

mpfr_t *_psi_tab;
mpfr_t *_sumlog_table;
mpfr_t _bernoulli[9];

mpfr_t _last_psi;
mpfr_t _psi_u;
mpfr_t _S_1;
mpfr_t _S_2;
mpfr_t _S_3;
mpfr_t _S_4;
mpfr_t _lambda;

void compute_small_psi(mpfr_t psi_value, long x);


void init_mpfr_vars() {
  mpfr_init_set_si(_tmp1,0, MPFR_RNDN);
  mpfr_init_set_si(_big_sumlog,0, MPFR_RNDN);
  mpfr_init_set_si(_small_sumlog,0, MPFR_RNDN);
  mpfr_init_set_si(_sumlog,0, MPFR_RNDN);
  mpfr_init_set_si(_mcarre,0, MPFR_RNDN);
  mpfr_init_set_si(_ncarre,0, MPFR_RNDN);
  mpfr_init_set_si(_n,0, MPFR_RNDN);
  mpfr_init_set_si(_m,0, MPFR_RNDN);
  mpfr_init_set_si(_n1,0, MPFR_RNDN);
  mpfr_init_set_si(_m1,0, MPFR_RNDN);
  mpfr_init_set_si(_logn,0, MPFR_RNDN);
  mpfr_init_set_si(_logm,0, MPFR_RNDN);
  mpfr_init_set_si(_den,0, MPFR_RNDN);
  mpfr_init_set_si(_ti,0, MPFR_RNDN);
  mpfr_init_set_si(_last_psi,0, MPFR_RNDN);
  mpfr_init_set_si(_psi_u,0, MPFR_RNDN);
  mpfr_init_set_si(_S_1,0, MPFR_RNDN);
  mpfr_init_set_si(_S_2,0, MPFR_RNDN);
  mpfr_init_set_si(_S_3,0, MPFR_RNDN);
  mpfr_init_set_si(_S_4,0, MPFR_RNDN);
  mpfr_init_set_si(_lambda,0, MPFR_RNDN);
}

void clear_mpfr_vars() {
  mpfr_clears(_tmp1, _big_sumlog, _small_sumlog, _sumlog, _mcarre, _ncarre, _n, _m, _n1, _m1, _logn, _logm, _den, _ti, _last_psi,\
	      _psi_u, _S_1, _S_2, _S_3, _S_4, _lambda, (mpfr_ptr) 0);
}

//inline LONG_DOUBLE LOG(long x) {return log(LONG_DOUBLE(x));}


template <class T_INT> inline T_INT MAX(T_INT a,T_INT b) {return (a>b)? a:b;}
template <class T_INT> inline T_INT MIN(T_INT a,T_INT b) {return (a<b)? a:b;}
template <class T_INT> inline T_INT ABS(T_INT a) {return (a>=0)? a:(-a);}
template <class T_INT> inline T_INT ODD(T_INT a) {return long(a)&1L;}


template <class T_INT> T_INT power(T_INT x,int n)
{
  T_INT r;
  if (ODD(n)) r=x;
  else        r=1;
  while ((n>>=1))
    {
      x *= x;
      if (ODD(n)) r*=x;
    }
  return r;
}

template <class T_INT> inline T_INT root(T_INT x,int n)
{
  T_INT u,v;
  u = 1;
  v = x/n+1;
  while (ABS(v-u)>1)
    {
      u=v;
      v=x;
      for (int i=0;i<n-1;i++) v/=u;
      v = (u*(n-1)+v)/n;
    }
  return MIN(u,v);
}

inline long root(long x,int n)
{
  long u,v;
  u = 1;
  v = (x-1)/n+1;
  while (ABS(v-u)>1)
    {
      u=v;
      v=x;
      for (int i=0;i<n-1;i++) v/=u;
      v = (u*(n-1)+v)/n;
    }
  return MIN(u,v);
}

int *_pi;
long pimax = 0;
void init_pi(long x)
{
  if (x<=0) {cerr << "bad x in init_pi, x = " << x << endl;exit(1);}
  _pi = new int[x+1];
  if (_pi==NULL) {cerr << "no memory for _pi !" << endl;exit(1);}
  long i;
  _pi[0] = 0; _pi[1] = 0; 
  for (i=2;i<=x;i++) _pi[i] = 1;
  long maxp = root(x,2);
  for (long p=2;p<=maxp;p++)
    {
      if (_pi[p])
	for (long m=p*p;m<=x;m+=p) _pi[m] = 0;
    }
  for (i=1;i<=x;i++) _pi[i] += _pi[i-1];
  pimax = _pi[x];
}

int *_p;
long pmax = 0;
long maxp = 0;
void init_p(long x)
{
  int i,pi;
  _p = new int[pimax+1];
  if (_p==NULL) {cerr << "no memory for _p !" << endl; exit(1);}
  _p[0] = 0;
  pi = 0;
  for (i=2;i<=x;i++) 
    if (_pi[i]-_pi[i-1]) _p[++pi] = i;
  pmax = _p[pimax];
  maxp = x;
}

long *_mu;
void init_mu(long x)
{
  _mu = new long[x+1];
  if (_mu==NULL) {cerr << "no memory for _mu !" << endl; exit(1);}
}

void build_mu(long a,long b,long *t=_mu)
{
  long i;
  long l = b-a;
  long sqrt_b = root(b,2);
  if (maxp<sqrt_b)
    {
      cerr << "Insufficient size of primetable" << endl;
      cerr << "maxp = " << maxp << ", sqrt(b) = " << sqrt_b << endl;
      exit(1);
    }
  for (i=0;i<l;i++) t[i] = 1;
  for (i=1;(i<=pimax)&&(_p[i]<=sqrt_b);i++)
    {
      long p = _p[i];
      long p2 = p; p2 *= p;
      // eliminate square factors
      long n  = -(a%p2);
      for (n=(n<0)?n+p2:n;n<l;n+=p2) t[long(n)] = 0;
      // compute t[]
      n = -(a%p);
      for (n=(n<0)?n+p:n;n<l;n+=p) t[long(n)] *= -p;
    }
  for (i=0;i<l;i++)
    {
      if (t[i]!=0)
	{
	  if (ABS(t[i])<(a+i)) t[i] = -t[i];
	  if (t[i]>0) t[i] =  1;
	  else        t[i] = -1;
	}
    }
}


void init_psi_tab(long l)
{
  _psi_tab = new mpfr_t[l+1];
  if (_psi_tab==NULL) {cerr << "no memory for _psi_tab !" << endl; exit(1);}
  for (int i=0; i <= l; i++) {
    mpfr_init_set_si(_psi_tab[i],i, MPFR_RNDN);
  }
}

void clear_psi_tab(long l) {
  for (int i=0; i <= l; i++)
    mpfr_clear(_psi_tab[i]);
  delete [] _psi_tab;
}


void build_psi(long a,long b)
{
  long n;
  long l = b-a;
  long sqrt_b = root(b,2);
  if (maxp<sqrt_b)
    {
      cerr << "Insufficient size of primetable" << endl;
      cerr << "maxp = " << maxp << ", sqrt(b) = " << sqrt_b << endl;
      exit(1);
    }
  for (n=0;n<l;n++) mpfr_set_si(_psi_tab[n],1, MPFR_RNDN);
  if (a==0) {
    mpfr_set_si(_psi_tab[0],0,MPFR_RNDN);
    mpfr_set_si(_psi_tab[1],0,MPFR_RNDN);
      }
  if (a==1)
    mpfr_set_si(_psi_tab[0],0,MPFR_RNDN);

  for (long k=1;(k<=pimax)&&(_p[k]<=sqrt_b);k++)
    {
      long p = _p[k];
      long m = MAX(p*p,a-(a%p));
      for (m = (m<a)?m+p:m; m<b; m+=p)	mpfr_set_si(_psi_tab[long(m-a)],0,MPFR_RNDN);
      long i=2;
      for (m=p*p;m<a;m*=p) i++;
      while (m<b)
	{
	  mpfr_set_si(_psi_tab[long(m-a)],i, MPFR_RNDN);
	  m*=p;i++;
	}
    }
  for (n=0;n<l;n++)
    {
      if (mpfr_zero_p(_psi_tab[n])) {
	mpfr_set(_psi_tab[n],_last_psi,MPFR_RNDN);
      }
      else  {
	// _psi_tab[n] = LOG(a+n)/_psi_tab[n]+lastpsi;
	mpfr_set_si(_tmp1, a+n, MPFR_RNDN);
	mpfr_log(_tmp1, _tmp1, MPFR_RNDN);
	mpfr_div(_tmp1,_tmp1,_psi_tab[n], MPFR_RNDN);
	mpfr_add(_psi_tab[n], _tmp1, _last_psi, MPFR_RNDN);
      }
      mpfr_set(_last_psi, _psi_tab[n], MPFR_RNDN);
    }
}

void init_bernoulli() {
  mpfr_init_set_si(_bernoulli[0],0,MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[1],1,MPFR_RNDN);
  mpfr_div_si(_bernoulli[1], _bernoulli[1], 6, MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[2],1,MPFR_RNDN);
  mpfr_div_si(_bernoulli[2], _bernoulli[2], 30, MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[3],1,MPFR_RNDN);
  mpfr_div_si(_bernoulli[3], _bernoulli[3], 42, MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[4],1,MPFR_RNDN);
  mpfr_div_si(_bernoulli[4], _bernoulli[4], 30, MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[5],5,MPFR_RNDN);
  mpfr_div_si(_bernoulli[5], _bernoulli[5], 66, MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[6],691,MPFR_RNDN);
  mpfr_div_si(_bernoulli[6], _bernoulli[6], 2730, MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[7],7,MPFR_RNDN);
  mpfr_div_si(_bernoulli[7], _bernoulli[7], 6, MPFR_RNDN);
  mpfr_init_set_si(_bernoulli[8],3617,MPFR_RNDN);
  mpfr_div_si(_bernoulli[8],_bernoulli[8], 510,MPFR_RNDN);
}

void clear_bernoulli() {
  for (int i=0; i < 9; i++)
    mpfr_clear(_bernoulli[i]);
}

long max_sumlog_table = 0;
void init_sumlog(long x)
{
  max_sumlog_table = x;
  _sumlog_table = new mpfr_t[x+1];
  mpfr_init_set_si(_sumlog_table[0],0,MPFR_RNDN);
  for (long i=1;i<=x;i++) {
    mpfr_init_set_si(_sumlog_table[i], i, MPFR_RNDN);
    mpfr_log(_sumlog_table[i], _sumlog_table[i], MPFR_RNDN);
    mpfr_add(_sumlog_table[i], _sumlog_table[i], _sumlog_table[i-1], MPFR_RNDN);
  }
}

void clear_sumlog(long x) {
  for (int i=0; i <= x; i++)
    mpfr_clear(_sumlog_table[i]);
  delete [] _sumlog_table;
}


//LONG_DOUBLE big_sumlog(long m, long n)
void big_sumlog(long m, long n)
{
  mpfr_set_si(_big_sumlog,0, MPFR_RNDN);
  mpfr_set_si(_n, n, MPFR_RNDN);
  mpfr_set_si(_m, m, MPFR_RNDN);
  mpfr_set_si(_ncarre, n, MPFR_RNDN); 
  mpfr_set_si(_mcarre, m, MPFR_RNDN);
  mpfr_mul_si(_ncarre, _ncarre, n, MPFR_RNDN);
  mpfr_mul_si(_mcarre, _mcarre, m, MPFR_RNDN);
  mpfr_set_si(_tmp1,1,MPFR_RNDN);
  mpfr_div(_mcarre,_tmp1, _mcarre, MPFR_RNDN);
  mpfr_div(_ncarre,_tmp1, _ncarre, MPFR_RNDN);

  mpfr_log(_logn, _n, MPFR_RNDN);
  mpfr_mul_d(_logn, _logn, n+0.5, MPFR_RNDN); 
  mpfr_log(_logm, _m, MPFR_RNDN);
  mpfr_mul_d(_logm, _logm, m-0.5, MPFR_RNDN);
  mpfr_add(_big_sumlog, _big_sumlog, _logn, MPFR_RNDN);
  mpfr_sub(_big_sumlog, _big_sumlog, _logm, MPFR_RNDN);
  mpfr_sub_si(_big_sumlog,_big_sumlog, n-m, MPFR_RNDN);

  long signe = -1 ;
  long imax  = 8;
  mpfr_set(_n1, _n, MPFR_RNDN);
  mpfr_set(_m1, _m, MPFR_RNDN);
  
  for (long i=1;i<=imax;i++)
    { 
      mpfr_mul(_n1, _n1, _ncarre, MPFR_RNDN);
      mpfr_mul(_m1, _m1, _mcarre, MPFR_RNDN);
      signe = -signe;
      mpfr_set_d(_den, signe * (2.0*i)*(i+i-1), MPFR_RNDN);
      mpfr_div(_ti,_bernoulli[i],_den, MPFR_RNDN);
      mpfr_sub(_tmp1, _n1, _m1, MPFR_RNDN);
      mpfr_mul(_ti, _ti, _tmp1, MPFR_RNDN);
      mpfr_add(_big_sumlog, _big_sumlog, _ti, MPFR_RNDN);
    }
}

//LONG_DOUBLE small_sumlog(long m, long n)
void small_sumlog(long m, long n)  
{
  mpfr_set(_small_sumlog, _sumlog_table[n], MPFR_RNDN);
  mpfr_sub(_small_sumlog, _small_sumlog, _sumlog_table[m-1], MPFR_RNDN);
}


//LONG_DOUBLE sumlog(long m, long n)
void sumlog(long m, long n)
{
  if (m>=max_sumlog_table) {
    mpfr_set(_sumlog, _big_sumlog, MPFR_RNDN);
    return;
  }
  
  if (n<=max_sumlog_table) {
    mpfr_set(_sumlog, _small_sumlog, MPFR_RNDN);
    return;
  }
  big_sumlog(max_sumlog_table+1,n);
  small_sumlog(m,max_sumlog_table);
  mpfr_add(_sumlog,_big_sumlog,_small_sumlog, MPFR_RNDN);
}

long *_xm;
long *_sqrt_xm;
void init_xm(long x,long u)
{
  _xm = new long[u+1];
  if (_xm==NULL) {cerr << "no memory for _xm !" << endl; exit(1);}
  _sqrt_xm = new long[u+1];
  if (_sqrt_xm==NULL) {cerr << "no memory for _sqrt_xm !" << endl; exit(1);}
  for (long m=1;m<=u;m++) {_xm[m] = x/m; _sqrt_xm[m]=root(_xm[m],2);}
}


void init(long x,long l,long u)
{
  init_pi(root(x/u,2));
  init_p(root(x/u,2));
  delete []_pi;
  init_mu(u+1);
  build_mu(0,u+1);
  init_xm(x,u);
  init_psi_tab(l);
  init_sumlog(100);
}

void test(long x,long l,long u)
{
  int i;
  for (i=0;i<root(x/u,2);i++) cout << "pi[" << i << "] = " << _pi[i]<< endl;
  for (i=0;i<root(x/u,2);i++) cout << "p[" << i << "] = " << _p[i]<< endl;
  for (i=0;i<=u;i++)          cout << "mu[" << i << "] = " << _mu[i] << endl;
  for (i=0;i<l;i++) {           
    cout << "psi[" << i << "] = ";
    mpfr_printf ("%.30Rf\n", _psi_tab[i]);
  }
}

long a = 0;
long b = 0;
inline long     mu(long n)     {return _mu[long(n)];}

void compute_lambda(long n) {
  mpfr_set(_lambda, _psi_tab[n-a], MPFR_RNDN);
  mpfr_sub(_lambda, _lambda, _psi_tab[n-a-1], MPFR_RNDN);
}

void compute_psi(mpfr_t psi_value, long x)
{
  if (x <= 800) {
    compute_small_psi(psi_value,x);
    return;
    }
  long  x3 = root(x,3);
  long   u = x3;
  long   l = u*8;
  init_mpfr_vars();
  mpfr_set_si(psi_value, 0, MPFR_RNDN);
  init_bernoulli();
  init(x,l,u);
  a = 1;
  b = a+l;
  mpfr_set_si(_last_psi,0,MPFR_RNDN);
  build_psi(a,b);
  mpfr_set(_last_psi, _psi_tab[long(b-1-a)], MPFR_RNDN);
  mpfr_set(_psi_u, _psi_tab[long(u-a)], MPFR_RNDN);
  mpfr_set(_S_1, _psi_u, MPFR_RNDN);
  //test(x,l,u);
  for (long n=1;n<=u;n++)
    {
      if (mu(n)==1)
	{
	  long xn = x/n;
	  //S_2 += sumlog(1,xn);
	  mpfr_set_si(_sumlog, 0, MPFR_RNDN);
	  sumlog(1,xn);
	  mpfr_add(_S_2, _S_2, _sumlog, MPFR_RNDN);
	  //_S_4 -= _psi_u*double(xn/u-u/n);
	  mpfr_mul_d(_tmp1, _psi_u, double(xn/u-u/n), MPFR_RNDN);
	  mpfr_sub(_S_4, _S_4, _tmp1, MPFR_RNDN);
	}
      if (mu(n)==-1)
	{
	  long xn = x/n;
	  //_S_2 -= sumlog(1,xn);
	  sumlog(1,xn);
	  mpfr_sub(_S_2, _S_2, _sumlog, MPFR_RNDN);
	  // _S_4 += _psi_u*double(xn/u-u/n);
	  mpfr_mul_d(_tmp1, _psi_u, double(xn/u-u/n), MPFR_RNDN);
	  mpfr_add(_S_4, _S_4,_tmp1, MPFR_RNDN);
	}
    }
  for (long m=2;m<=u;m++)
    {
      compute_lambda(m);
      if (! mpfr_zero_p(_lambda))
	{
	  long xm = _xm[m];
	  long S  = 0;
	  for (long n=1;n<=u;n++)
	    {
	      if (mu(n)==1)  S += xm/n;
	      if (mu(n)==-1) S -= xm/n;
	    }
	  //_S_3 += lambda(m)*S;
	  mpfr_set_si(_tmp1, S, MPFR_RNDN);
	  mpfr_mul_si(_lambda, _lambda, S, MPFR_RNDN);
	  mpfr_add(_S_3, _S_3, _lambda, MPFR_RNDN);
	}
    }

  long xu = x/u;
  //cout << "u= " << u << "   xu= " << xu << endl;
  while (a<xu)
    {
      //cout << "a= " << a << endl;
      for (long m = 1; m <= u; m++)
	{
	  if (mu(m)==0) continue;
	  long xm = _xm[m];
	  long n1 = MAX(long(u/m),xm/b)+1;
	  long n = n1;
	  long n2;
	  n2 = MIN(xm/a,long(xm/u));
	  if (mu(m)==1) 
	    {
	      while (n <= n2)
		{
		  long q = xm/n;
		  long nextn = xm/q + 1;
		  if (nextn <= n2) {
		    // _S_4 += psi(q)*LONG_DOUBLE(nextn-n);
		    mpfr_set_si(_tmp1, nextn-n, MPFR_RNDN);
		    mpfr_mul(_tmp1, _tmp1, _psi_tab[q-a], MPFR_RNDN);
		    mpfr_add(_S_4, _S_4, _tmp1, MPFR_RNDN);
		  }
		  else {
		    //_S_4 += psi(q)*LONG_DOUBLE(n2-n+1);
		    mpfr_set_si(_tmp1, n2-n+1, MPFR_RNDN);
		    mpfr_mul(_tmp1, _tmp1, _psi_tab[q-a], MPFR_RNDN);
		    mpfr_add(_S_4, _S_4, _tmp1, MPFR_RNDN);
		  }
		  n = nextn;
		}
	    }
	  else
	    {
	      while (n <= n2)
		{
		  long q = xm/n;
		  long nextn = xm/q + 1;
		  if (nextn <= n2) {
		    // _S_4 -= psi(q)*LONG_DOUBLE(nextn-n);
		    mpfr_set_si(_tmp1, nextn-n, MPFR_RNDN);
		    mpfr_mul(_tmp1, _tmp1, _psi_tab[q-a], MPFR_RNDN);
		    mpfr_sub(_S_4, _S_4, _tmp1, MPFR_RNDN);
		  }
		  else {
		    // _S_4 -= psi(q)*LONG_DOUBLE(n2-n+1);
		    mpfr_set_si(_tmp1, n2-n+1, MPFR_RNDN);
		    mpfr_mul(_tmp1, _tmp1, _psi_tab[q-a], MPFR_RNDN);
		    mpfr_sub(_S_4, _S_4, _tmp1, MPFR_RNDN);
		  }
		  n = nextn;
		}
	    }
	}
      //mpfr_printf("_S_4= %.30Rf\n", _S_4);
      a+=l;
      if (a>=xu) break;
      b = MIN(long(a+l),xu);
      build_psi(a,b);
      //last_psi = psi(b-1);
      mpfr_set(_last_psi, _psi_tab[b-1-a], MPFR_RNDN);
    }
  mpfr_add(psi_value, psi_value, _S_1, MPFR_RNDN);
  mpfr_add(psi_value, psi_value, _S_2, MPFR_RNDN);
  mpfr_sub(psi_value, psi_value, _S_3, MPFR_RNDN);
  mpfr_sub(psi_value, psi_value, _S_4, MPFR_RNDN);
  clear_mpfr_vars();
  clear_psi_tab(l);
  clear_sumlog(100);
  clear_bernoulli();
  delete[] _mu;
  delete[] _p;
  delete[] _xm;
  delete[] _sqrt_xm;
}

const int small_primes[139]={2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139,
149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223,
227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293,
307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383,
389, 397, 401, 409, 419, 421, 431, 433, 439, 443, 449, 457, 461, 463,
467, 479, 487, 491, 499, 503, 509, 521, 523, 541, 547, 557, 563, 569,
571, 577, 587, 593, 599, 601, 607, 613, 617, 619, 631, 641, 643, 647,
653, 659, 661, 673, 677, 683, 691, 701, 709, 719, 727, 733, 739, 743,
751, 757, 761, 769, 773, 787, 797};

void compute_small_psi(mpfr_t psi_value, long x) {
  //cout << "In compute small psi x= " << x << endl;
  mpfr_t pmpfr, logp;
  mpfr_inits2(128,pmpfr,logp, (mpfr_ptr) 0);
  mpfr_set_si(psi_value,0,MPFR_RNDN);
  int i=0;
  long p=small_primes[i++];
  while (p <= x) {
    mpfr_set_si(pmpfr,p, MPFR_RNDN);
    mpfr_log(logp, pmpfr, MPFR_RNDN);
    long q=p;
    while (q <= x) {
      mpfr_add(psi_value, psi_value, logp, MPFR_RNDN);
      q *= p;
    }
    p=small_primes[i++];
  }
  mpfr_clears(pmpfr, logp, (mpfr_ptr) 0);
}
