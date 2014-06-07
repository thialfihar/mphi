#ifndef __FUNCTIONS_H__
#define __FUNCTIONS_H__
#include <cstdio>
#include <cstring>

#include <algorithm>
#include <vector>
#include <gmpxx.h>

using std::pair;
using std::vector;

typedef struct {
    mpz_class prime;
    unsigned int exponent;
} Factor;
typedef vector<Factor> Factors;

Factors factorize(mpz_class n);
vector<mpz_class> get_divisors(mpz_class n);
mpz_class phi_n(unsigned int n, unsigned int k);
unsigned int gcd(unsigned int a, unsigned int b);

#endif
