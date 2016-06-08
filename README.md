# Computing h(n), the maximum value of a product
# of distinct primes whose sum is not greater than n.

For example, h(12) = 42 because  the maximum product
of _distinct primes_ whose sum is not larger than 12 is  2 x 3 x 7 = 42.


The h function is rapidly increasing, h(n) is an integer whose
number of digits is about 0.43 sqrt(n log n), and it is impossible
to write the decimal expansion of h(n) for moderate values of n.

Nethertheless it is possible to describe the factorization of h(n).

### The factorization of h(n)

Let n be a positive integer. We define p<sub>k</sub> as
the largest prime p such that s<sub>k</sub>= 2 + 3 + .... + p<sub>k</sub> <= n,
and we denote dk  the difference n - s<sub>k</sub>.

For example, for n=12 and N=10^9, we have, respectively,

| n  | p<sub>k</sub> |  s<sub>k</sub>  | d<sub>k</sub> |
| ---------: | ---------: | ---------: | --------: |
| 12  | 5 | 10  | 2 |
| 10^9 | 151,057  | 24,739,512,092,254,535 | 133,215 |

Let us define N<sub>k</sub> = 2 x 3 x .... x p<sub>k</sub>. The value
h(n) **is rather close to N<sub>k</sub>**.

More precisely, h(n) may be written as a product
h(n) = N<sub>k</sub>  G(p<sub>k</sub>, d<sub>k</sub>)
where  G(p,m) is a rational number which is defined in
[Landau's function for one million billions] (https://eudml.org/doc/10854).
This rational number  G(p,m)  has the following property:

There exists a small integer s and primes q<sub>j</sub> and Q<sub>j</sub> such that
- q<sub>1</sub> q<sub>2</sub> < ... < q<sub>s</sub> <= p<sub>k</sub>
   < Q<sub>1</sub> < Q<sub>2</sub> < ... < Q<sub>s</sub>
- G(p<sub>k</sub>, d<sub>k</sub>) = (Q<sub>1</sub>  Q<sub>2</sub>  ...
  Q<sub>s</sub>) / (q<sub>1</sub>  q<sub>2</sub>  ...  q<sub>s</sub>).

More over, all the Q<sub>j</sub>, and all of the q<sub>j</sub>, 
__except perhaps q<sub>1</sub>_, are close to p<sub>k</sub>.  

For n=12 and N=10^9, we have, respectively,

|    n  | p<sub>k</sub> |  q<sub>1</sub>... q<sub>s</sub> | Q<sub>1</sub>... Q<sub>s</sub> |G( p<sub>k</sub> , d<sub>k</sub> ) |
| ---------: | ---------: | ---------: | ---------: | :--------: |
|12  |  5 | 5 | 7 | 7/5|
|10<sup>9</sup> | 151,057 | 17,881 | 151,091|  151,091 / 17,881|


- For n=12, N<sub>k</sub>=2 x 3 + 5 = 30,  G(5,2) = 7/5 and h(n) = 30 x 7 /5 = 42.

- For n=10^9, N<sub>k</sub> is the product of all the primes uo to 151,057. Its decimal
expansion has 65,449 digits, while the expansion of h(n) has 65,450 digits.
But the fraction h(n)/N<sub>k</sub>=G(151,057, 133,215) is reduced to 151,091/17,881.


# Command-line  hcompute
-------------------------------
```
Usage: hcompute  [OPTION] n
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


# Prerequisites

Before using this repository you need to install a  _primesum_
function on your computer. This function computes the sum of
primes upto x, doing  O(x^(2/3)/log^2 x) elemantary operations.

You will find  a very good **primesum**  command on the repository
[primesum](https://github.com/kimwalisch/primesum)
of Kim Walisch. This implementation uses threads and this may easily
divide by 10 or more the time of computation on a computer with a large
nymber of cores.
The command primesum must be in the directory where you download Gh, or
must be in your $PATH.

Build instructions (Unix-like OSes)
-----------------------------------
You need to have installed a C++ compiler and GNU make.

Download
[Gh.tar.gz](https://dl.bintray.com/mhdeleglise/Gh/Gh.tar.gz)
and build it using:

```
$ ./make
```

References
----------
1. M. Deleglise, J. Rivat,
[ Computing pi(x): The Meissel, Lehmer,Lagarias, Miller, Odlyzko method]
(http://www.ams.org/journals/mcom/1996-65-213/S0025-5718-96-00674-6/)
. Math. Comp., 65 (1996), no. 213, 235-245.
2. M. Deleglise and J.-L. Nicolas, [Maximal product of primes whose sum
is bounded](http://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=spm&paperid=45&option_lang=eng),
Bull. Proc. of the Steklov Institute 17 (2013), 82-112)

3. M. Delegise, J.-L. Nicolas and P. Zimmermann,
[Landau's function  for one million billions](https://eudml.org/doc/10854)
, J. Theor. Nombres Bordeaux. 20 (2008), no. 3, 625-671.
