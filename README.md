# Computing h(n)

##### the maximum value of a product of distinct primes whose sum is not greater than n.

The function h(n) is defined in 
[Maximal Product of Primes Whose Sum is Bounded]
(http://www.mathnet.ru/php/archive.phtml?wshow=paper&jrnid=spm&paperid=45&option_lang=eng)


For example, h(12) = 42 because  the maximum value of a product
of primes whose sum is not larger than 12 is  2 x 3 x 7 = 42.

Let us denote p(j) the  j em prime,

- s(j) =  2 + 3 + .... + p(j)  the sum of the j first primes

- N(j)=   2 x 3 + .... x p(j)  the product of the j first primes

Let n be a positive integer. We define k(n),
the integer k such that s(k) <= n < s(k+1)

and we denote d(k) the difference d(k) = n - s(k).

The value h(n) may be written as a product h(n) = Nk x G(pk, dk)

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
