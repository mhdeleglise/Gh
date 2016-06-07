# Computing h(n)

##### the maximum value of a product of distinct primes whose sum is not greater than n.

The function h(n) is defined in 
[Maximal Product of Primes Whose Sum is Bounded]
(http://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=spm&paperid=45&option_lang=eng).
For example, h(12) = 42 because  the maximum value of a product
of _distinct primes_ whose sum is not larger than 12 is  2 x 3 x 7 = 42.

Let n be a positive integer. We define pk as
the largest prime p such that sk= 2 + 3 + .... + pk <= n,
and we denote dk  the difference n - sk.

For example, for n=12 and N=10^9, we have, respectively,

| n  | pk |  sk  | dk |
| ---------: | ---------: | ---------: | --------: |
| 12  | 5 | 10  | 2 |
| 10^9 | 151057  | 24739512092254535 | 133215 |

Let us define Nk = 2 x 3 x .... x pk. The value h(n) is relatively
close to Nk. More precisely, h(n) may be written as a product h(n) = Nk x G(pk, dk)
where  G(p,m) is a rational number which is defined in
[Landau's function for one million billions] (https://eudml.org/doc/10854).
This rational number  G(p,m)  has the following property:

there exists a small integer s and primes qj and Qj such that
q1 < q2 < ... < qs <= pk < Q1 < Q2 < ... Qs
- Num = Q1 x Q2  ...  x Qs
- Den  =  q1 x q2  ...  x qs
and

- G(pk,dk) = Num/Den.


h(n) is an huge number whose number of digits is about
0.43 x sqrt(n log n), and it is impossible to write the decimal
expansion of h(n).

Let us come back Coming  to the example n=1000000000. In this case the
number Nk is the product of all the primes uo to 151057. Its decimal
expansion has 65450 digits. But the fraction G(pk, dk) = G(151057,
133215) is reduced to 151091/17881.




```
hcompute 100000000000000000000
```

```
hcompute -v 100000000000000000000
```

```
hcompute -p 100000000000000000000
```

```
hcompute -l 100000000000000000000
```
