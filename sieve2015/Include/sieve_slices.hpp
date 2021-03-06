// -*- C++ -*-
#ifndef sieve_by_slice_hpp
#define sieve_by_slice_hpp
#include"mylib.h"
#include"primes.h"
#include"prime_generator.h"

template<class btable, class longint> void 
sieve_by_slice<btable, longint>::create(int k0, int l0, long64 slice_size, 
					longint startx, int presieve_base, sieve_type t)
{
  sieve_t= t;
  k = k0; l = l0;
  access_frame::create(presieve_base, k, l);
  // Attention
  window_size=next_mult<longint>(slice_size, get_period());
  window_start=(startx/window_size)*window_size;
  window_end=window_start+window_size-1;
  /*
  window_start=prev_mult<longint>(startx, get_period());
  window_end=next_mult<longint>(startx+slice_size, get_period())-1;
  window_size=(window_end-window_start+1);
  */
  table_inverses.create((int)get_period());
  //cout << "window_size= " << window_size << endl;
  //cout << "lower_index64 de window_size= " << lower_index64(window_size) << endl;
  btable::create(1 + lower_index64(window_size));
  btable::unset_bit(0);
  last_total = 0;
  if (sieve_t == AUTO_SIEVE)
    {
#ifdef DEBUG_SV
      cout << "sieve_by_slice<btable>::create: AUTO_SIEVE is set: Eratosthenes_sieve will be done\n";
#endif
      eratosthenes();
    }
  //display();
}

template<class btable, class longint> void 
sieve_by_slice<btable, longint>::set_around(longint x)
{
  //cout << "In set_around x= " << x << " wstart= " << window_start << "    wend= " << window_end << "  window_size= " << window_size <<  endl;
  if ((x < window_start) || (x > window_end)) {
    long q = x / window_size;
    window_start = q * window_size;
    window_end   = window_start + window_size-1;
    btable::fill();
    btable::unset_bit(0);
    if (sieve_t == AUTO_SIEVE)
      {
#ifdef DEBUG_SV
	cout << "sieve_by_slice<btable>::create: AUTO_SIEVE is set: Eratosthenes_sieve will be done\n";
#endif
	eratosthenes();
      }
    //cout << "Now  x = " << x << "   window_start= " << window_start << "    end= " << window_end << endl;
  }
}


template<class btable, class longint> void 
sieve_by_slice<btable, longint>::sieve_by(long64 p)
{
  long64 inc = p<<1;
  long64 true_inc = (p<<2) + (p<<1);
  
  int r  = p % 6;
  longint start;
  long64 index;
  
  if (r==1)
    {
      // j=0
      start = p;
      if (start < window_start)
	start += ((window_start-start)/true_inc + 1) * true_inc;
      index = lower_index64(start-window_start);
      for ( ; index < btable::get_bit_size(); index += inc)
	btable::unset_bit(index);

      // j=1
      start= p+(p<<2);   // 5p
      if (start < window_start)
	start += ((window_start-start)/true_inc + 1) * true_inc;
      index = lower_index64(start-window_start);
      for ( ; index < btable::get_bit_size(); index += inc)
	btable::unset_bit(index);
    }
  else // r=5
    {
      //j=0;
      start= p+(p<<2);
      if (start < window_start)
	start += ((window_start-start)/true_inc + 1) * true_inc;
      index = lower_index64(start-window_start);
      for ( ; index < btable::get_bit_size(); index += inc)
	btable::unset_bit(index);

      //j=1
      start= p;
      if (start < window_start)
	start += ((window_start-start)/true_inc + 1) * true_inc;
      index = lower_index64(start-window_start);
      for ( ; index < btable::get_bit_size(); index += inc)
	btable::unset_bit(index);
    }
}


template<class btable, class longint> void 
sieve_by_slice<btable, longint>::shift_window_forward()
{

#ifdef DEBUG_SV
  cout << "In sieve_slice::shift_window_forward forward start = " << window_start << endl;
  cout << "                                           end  = " << window_start+window_size-1 << endl;
#endif
  
  if (btable::init_counters()) {
    if (sieve_t==AUTO_SIEVE)
      last_total = count(window_start+window_size-1);
  }
  window_start += window_size;
  window_end += window_size;
  btable::fill();
  btable::unset_bit(0);
  if (sieve_t == AUTO_SIEVE)
    {
#ifdef DEBUG_SV
      cout << "sieve_slice::forward, AUTO_SIEVE is set: Eratosthenes sieve will be done" << endl;
#endif
      eratosthenes();
    }
}

template<class btable, class longint> void 
sieve_by_slice<btable, longint>::shift_window_backward()
{
  window_start -= window_size;
  window_end   -= window_size;
  btable::fill();
  btable::unset_bit(0);
  if (sieve_t == AUTO_SIEVE)
    {
      eratosthenes();
      last_total -= btable::count((int)lower_index64(window_size-1));  
    }
}

template<class btable, class longint> void 
sieve_by_slice<btable, longint>::display(int how_many)
{
  cout << "Sieve_linear::display size= " << window_size << "   [" << window_start << ", " << window_start + window_size - 1 << "]"  \
       << "    sieve_type= " << sieve_t;
  cout << "    last_total= " << last_total << endl;
  cout << "    btable::get_bit_size= " << btable::get_bit_size() << endl;
  int cnte=0;
  for (int i = 1; i < btable::get_bit_size(); i++)
    if (btable::get_bit(i))
      {
	cnte++;
	cout.width(8);
	cout << get_integer(i);
	if (cnte==how_many)
	  break;
      }
  cout << "\n\n";
}

template<class btable, class longint> void 
sieve_by_slice<btable, longint>::display_counts() {
    cout << "Sieve_linear::display(): Les comptes\n";
    for ( int i = 0 ; i < btable::bit_size ; i++)
      {
	cout.width(8);
	cout << count(get_integer(i));
      }
    cout << "\n\n";
}

template<class btable, class longint> void 
sieve_by_slice<btable, longint>::eratosthenes()
{
#ifdef DEBUG_SV
  cout << "In ERATOS\n";
#endif
  if (!presieved_primes::primes_initialized) {
    cout << "Error in sieve_by_slice::eratosthenes:   prime_table is not initailized\n";
    error();
  }
  if (sieve_t == NO_SIEVE)
    {
      cout << "Eratosthenes is not defined for this bit table!\n";
      error();
    }
  long64 p;
  long64 ip =1;

  if (!window_start && (l == 1)) btable::unset_bit(1);
  // bound for the primes we need
  double bound_p = sqrt((double)get_window_end()+2);
  // primes of primes:: cannot exceed max_prime()
  double bound1_p = min(double(presieved_primes::max_prime()),bound_p);
  while (p = presieved_primes::prime(ip++), p <= bound1_p)
   {
     sieve_by(p);
     if (p >= window_start && mod32(p, k) == l) 
       btable::set_bit(lower_index64(p-window_start));
   }
  if (bound_p > bound1_p) {
    cerr << "presieved_primes::max_prime() = " << presieved_primes::max_prime()\
	 << "  is too small to sieve until window_end = " << get_window_end() << endl;
    cerr << "We finish eratosthenes_sieve using a prime_generator " << endl;
    prime_generator pg(long(min(1.1*bound_p, 4000000000.0)), presieved_primes::max_prime());
    //pg.display();
    p=pg.next_prime();
    while (p < bound_p) {
      sieve_by(p);
      if (p >= window_start) 
	btable::set_bit(lower_index64(p-window_start)); 
      p=pg.next_prime();
    }
  }
  // cout << "\nTous les sieves sont faits. Faut il réinitialiser des compteurs\n";
  btable::init_counters();
  left_index = 0;
  right_index = btable::get_bit_size();
#ifdef DEBUG_SV
  cout << "OUT ERATOS\n";
#endif
}

template<class btable, class longint> longint
sieve_by_slice<btable, longint>::get_first_prime()
{
  left_index=0;
  for (;;)
    {
      while(++left_index < btable::get_bit_size())
	{
	  if (btable::get_bit(left_index))
	    return get_integer(left_index);
	}
      shift_window_forward();
    }
}

template<class btable, class longint> longint
sieve_by_slice<btable, longint>::get_next_prime()
{
  //cout << "In get_next_prime:  left_index = " << left_index << "    image " << get_integer(left_index) << endl;
  if (sieve_t == NO_SIEVE)
    {
      cout << "sieve_by_slice de type NO_SIEVE\n";
      cout << "sieve_by_slice::get_next_prime n'est pas défini\n";
      error();
    }
  for (;;)
    {
      while(++left_index < btable::get_bit_size())
	{
	  //cout << "left_index = " << left_index << "    image " << get_integer(left_index) << endl;
	  if (btable::get_bit(left_index))
	    return get_integer(left_index);
	  //else cout << "not prime " << endl;
	}
      //cout << "sieve_by_slices::forward\n";
      shift_window_forward();
    }
}

template<class btable, class longint> longint
sieve_by_slice<btable, longint>::get_next_prime_without_shifting()
{
  for (;;)
    {
      while(++left_index < btable::get_bit_size())
	{
	  if (btable::get_bit(left_index))
	    return get_integer(left_index);
	}
      return 0;
    }
}

template<class btable, class longint> longint  
sieve_by_slice<btable, longint>::get_previous_prime()
{
  //cout << "In get_previous_prime\n";
  //cout << " wstart= " << window_start << "    wend= " << window_end << "  window_size= " << window_size <<  endl;
  do
    {
      while (right_index--)
	{
	  if (btable::get_bit(right_index)) {
	    //cout << " wstart= " << window_start << "    wend= " << window_end << "  window_size= " << window_size <<  endl;
	    return get_integer(right_index);
	  }
	}
      if (window_start) 
	{
	  //cout << "backward\n";
	  shift_window_backward();
	}
    } while (right_index);
  cout << "ZZZZZZZZZZ\n";
  return 0;
}

template<class btable, class longint> void sieve_by_slice<btable, longint>::set_indexes(longint x)
{
  set_around(x);

  left_index = lower_index64(x - window_start);

  right_index = (int)lower_index64(x - window_start);
  if (get_integer(right_index) < x)
    right_index++;
}


template<class btable, class longint> longint
sieve_by_slice<btable, longint>::get_previous_prime(longint x)
{
  set_around(x);
  //cout << "After set_around x= " << x << " wstart= " << window_start << "    wend= " << window_end << "  window_size= " << window_size <<  endl;
  //display();
  right_index = (int)lower_index64(x - window_start);
  if (get_integer(right_index) < x)
    right_index++;
  //cout << "right_index= " << right_index << "  integer = " << get_integer(right_index) << endl;
  return get_previous_prime();
}

template<class btable, class longint> longint
sieve_by_slice<btable, longint>::get_next_prime(longint x)
{
  //set_around(x);
  //display();
  left_index = (int)lower_index64(x - window_start);
  //cout << "left_index= " << left_index << "  integer = " << get_integer(left_index) << endl;
  return get_next_prime();
}


template<class btable, class longint> longint
sieve_by_slice<btable, longint>::count(longint x)
{
  if (sieve_t == AUTO_SIEVE)
    while (x >= window_start + window_size) 
      {
	shift_window_forward();
      }

  return(last_total + 
	 btable::count((int)lower_index64(x-window_start))); 
}

template<class btable, class longint> int
sieve_by_slice<btable, longint>::belong_to_window(longint x)
{
  return (x >= window_start) && (x < window_start + window_size);
}

#endif
