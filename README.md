# Computing h(n)

## the maximum value of a product of distinct primes
## whose sum is not greater than n.

For example, h(10) = 2 x 3 x 5 = 30 is the maximum value of a product
of primes whose sum is not larger than 10.

Let us denote p(j) the  j em prime,

- s(j) =  2 + 3 + .... + p(j)  the sum of the j first primes

- N(j)=   2 x 3 + .... x p(j)  the product of the j first primes

Let n be a positive integer. We define k(n),
the integer k such that s(k) <= n < s(k+1)

and we denote d(k) the difference d(k) = n - s(k).

The value h(n) may be written as a product h(n) = Nk x G(pk, dk)

Some basic Git commands are:
```
hcompute 100000000000000000000
```

```
hcompute -v 100000000000000000000
```

```
hcompute -P 100000000000000000000
```


```
hcompute -l 100000000000000000000
```
