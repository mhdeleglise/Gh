# Computing h(n)

#### Maximum value of a product of distinct primes whose sum is not greater than n.

The function h(n) is defined in 
[Maximal Product of Primes Whose Sum is Bounded]
(http://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=spm&paperid=45&option_lang=eng).


For example, h(12) = 42 because  the maximum value of a product
of _distinct primes_ whose sum is not larger than 12 is  2 x 3 x 7 = 42.


The h function is rapidly increasing, h(n) is an integer whose
number of digits is about 0.43 x sqrt(n log n), and it is impossible
to write the decimal expansion of h(n) for moderate values of n.

Nethertheless it is possible to describe the factorization of h(n).

#### The factorization of h(n)

Let n be a positive integer. We define pk as
the largest prime p such that sk= 2 + 3 + .... + pk <= n,
and we denote dk  the difference n - sk.

For example, for n=12 and N=10^9, we have, respectively,

| n  | pk |  sk  | dk |
| ---------: | ---------: | ---------: | --------: |
| 12  | 5 | 10  | 2 |
| 10^9 | 151057  | 24739512092254535 | 133215 |

Let us define Nk = 2 x 3 x .... x pk. The value h(n) is relatively
close to Nk.

More precisely, h(n) may be written as a product h(n) = Nk x G(pk, dk)
where  G(p,m) is a rational number which is defined in
[Landau's function for one million billions] (https://eudml.org/doc/10854).
This rational number  G(p,m)  has the following property:

There exists a small integer s and primes qj and Qj such that
- q1 < q2 < ... < qs <= pk < Q1 < Q2 < ... Qs
- G(pk,dk) = (Q1 x Q2  ...  x Qs)/(Q1 x Q2  ...  x Qs).

| n  | q1 q2  ...  qk  | Q1 Q2 ... Qk |  G(pk,dk)  | 
| ---------: | ---------: | ---------: | :--------: |
| 12  | 5 | 7  | 7/5 |
| 10^9 | 151091  | 17881 | 151091/17881 |

For n=10^9, Nk is the product of all the primes uo to 151057. Its decimal
expansion has 65449 digits, while the expansion of h(n) has 65450 digits.
But the fraction h(n)/Nk=G(151057, 133215) is reduced to 151091/17881.

# Prerequisites

Before using this repository you need to install the _primesum_
function on your computer. This function computes the sum of
primes upto x, doing  O(x^(2/3)/log^2 x) elemantary operations.
The Kim Walisch implementation use threads, and, on a computer
with a large numbers of cores, this time may be easily divided by 10.

You will find  it on the repository
[primesum](https://github.com/kimwalisch/primesum)
of Kim Walisch.
And this function must be in the directory where you download Gh, on
must be in your $PATH.

# Command-line  hcompute
-------------------------------
```
Usage: primesum [OPTION] n
Computes the factorization of h(n) for n <= 10^35 and some related values.

Options:

  -f  prints the factorization of h(n)
  -l, prints the log of h(n) 
  -p, print the largest prime factor of h(n)
  -v, print h(n) if it has no more than 100000 digits.
```

# Command-line skcompute
----------------------

This repository gives you the skcompute function

```
skcompute n
```
which prints the values of pk and sk corresponding to n, for n upto 10^35.
