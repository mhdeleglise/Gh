# Computing h(n), the maximum value of a product
# of distinct primes whose sum is not greater than n.

For example, h(12) = 42 because  the maximum product
of _distinct primes_ whose sum is not larger than 12 is  2 x 3 x 7 = 42.


The h function is rapidly increasing, h(n) is an integer whose
number of digits is about 0.43 sqrt(n log n), and it is impossible
to write the decimal expansion of h(n) for moderate values of n.

However it is possible to describe the factorization of h(n).

### The factorization of h(n)

Let n be a positive integer. We define p<sub>k</sub> as
the largest prime p such that s<sub>k</sub>= 2 + 3 + .... + p<sub>k</sub> <= n,
and we denote d<sub>k</sub>  the difference n - s<sub>k</sub>.

For example, for n=12, 213, 10^9, we have, respectively,

| n  | p<sub>k</sub> |  s<sub>k</sub>  | d<sub>k</sub> |
| ---------: | ---------: | ---------: | --------: |
| 12   | 5    | 10    | 2 |
| 213 | 37  | 197  | 16 |
| 10^9 | 151,057  | 24,739,512,092,254,535 | 133,215 |

Let us define N<sub>k</sub> = 2 x 3 x .... x p<sub>k</sub>, that we
will write write more concisely [ 2-p<sub>k</sub> ]. The value
h(n) **is rather close to N<sub>k</sub>**.

More precisely, h(n) may be written as a product
h(n) = N<sub>k</sub>  G(p<sub>k</sub>, d<sub>k</sub>)
where  G(p,m) is a rational number which is defined in
[Landau's function for one million billions] (https://eudml.org/doc/10854).
This rational number  G(p,m)  has the following property:

There exists a small integer s and primes q<sub>j</sub> and Q<sub>j</sub> such that
- q<sub>1</sub> <  q<sub>2</sub> < ... < q<sub>s</sub> <= p<sub>k</sub>
   < Q<sub>1</sub> < Q<sub>2</sub> < ... < Q<sub>s</sub>
- G(p, m) = (Q<sub>1</sub>  Q<sub>2</sub>  ...
  Q<sub>s</sub>) / (q<sub>1</sub>  q<sub>2</sub>  ...  q<sub>s</sub>).

More over, all the Q<sub>j</sub>, and all of the q<sub>j</sub>, 
**except perhaps q<sub>1</sub>**, are close to p<sub>k</sub>.  

For n=12, 213, 10^9,  we have, respectively,

|    n  | p<sub>k</sub> |  d<sub>k</sub> | q<sub>1</sub>... q<sub>s</sub> | Q<sub>1</sub>... Q<sub>s</sub> |G( p<sub>k</sub> , d<sub>k</sub> )
| ---------:  | ---------: |  ---------: | ---------: | ---------: | :--------: |
|12   |  5  | 2   | 5 | 7 | 7/5|
|213 | 37 | 16 | 31  37 | 41  43 | (41 x 43) / (31 x 37)|
|10<sup>9</sup> | 151,057 | 133,251 |17,881 | 151,091|  151,091 / 17,881|

- For n=12, [2-5] =2 x 3 x 5,  p<sub>k</sub> = 5,  d<sub>k</sub> =
  2,  G(5,2) = 7/5 and h(n) = 2 x 3 x 5  x 7 /5  =42.

- For n=213, G(37, 16) = (41 x 43 ) / (37 x 31) , thus h(n) = [ 2-37 ]  x 41 x 43 / 37 / 31;

- For n=10^9, the
decimal expansion of N<sub>k</sub> has 65,449 digits, while the
expansion of h(n) has 65,450 digits. However h(n) is concisely
described by its factorization [ 2-151,057] x 151,091 / 17,881.



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

Caution: computing log h(n) is hard. The computation needs a few
seconds for n=10^20, and, each time you multiply n by 10, you
approximatively double the time of computation of log(h(n)).
The essential part of this
computation is to get the value theta(p<sub>k</sub>).
Theta function  is defined by theta(x) = sum{log p, for p prime, p<=x}.
The Meissel method doesn't work for this sum. In (4) an
algorithm computing theta(x) in time O(x^(2/3+epsilon) is given.
It is used in this package.

# Command-line skcompute
----------------------

This repository gives you the skcompute function, which prints out the
values p<sub>k</sub> an <sub>k</sub> associated to n.

```
skcompute n
```
which prints the values of pk and sk corresponding to n, for n upto 10^35.


# Prerequisites

Before using this repository you need to install a  _primesum_
function on your computer. This function computes the sum of
primes upto x, doing  O(x^(2/3)/log^2 x) elemantary operations.

You will find  a very efficient **primesum**  command on the repository
[primesum](https://github.com/kimwalisch/primesum)
of Kim Walisch. This implementation uses threads and this may easily
divide by 10 or more the time of computation on a computer with a large
number of cores.
The command primesum must be put  in the directory where you download Gh, or
in one of the directories figuring in the value of your PATH environment variable.

Build instructions (Unix-like OSes)
-----------------------------------
You need to have installed a C++ compiler and GNU make.
It works on PC linux machines and on my macbook pro.

Download this repository, and build it using

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

4. M. Deleglise, J. Rivat, [Computing Psi(x)](http://www.ams.org/journals/mcom/1998-67-224/S0025-5718-98-00977-6/).
Math. Comp. 67 (1998), no. 224, 1691-1696].

