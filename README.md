# Computing h(n)

## the maximum value of a product of distinct primes
## whose sum is not greater than n.

For example, h(10) = 2 x 3 x 5 = 30 is the maxumum value of a product
of primes whose sum is not larger than 10.

Let us denote p_j the  jeme prime,

- s_j =  2 + 3 + .... + p_j  the sum of the j first primes

- N_j=   2 x 3 + .... x p_j  the product of the j first primes

Let n be a positive integer. We define k(n),
the integer k such that

s_k <= n < s_(k+1)

and

we denote

d_k the diffeence d_k = n - s_k.

The value h(n) may be written as a product

h(n) = N_k x G(p_k, delta_k)

