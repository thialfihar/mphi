#!/usr/bin/env python
import math
import subprocess
from sympy import symbols, lcm

def factor(n):
    factors = {}
    for p in subprocess.check_output(["factor", str(n)]).split(':')[1].strip().split(' '):
        p = int(p)
        factors[p] = factors.get(p, 0) + 1

    return factors.items()

def phi_n(n, k):
    factors = factor(k)
    result = 1
    for i in range(n):
        for (p, e) in factors:
            result *= p ** (e * n) - p ** (e * n - i - 1)

    return result

def lambda_n(n, k):
    factors = factor(k)
    result = 1
    x = symbols('x')
    term = 1
    for i in range(1, n + 1):
        term = lcm(term, x ** i - 1)

    #print term
    for p, e in factors:
        result = lcm(result, int(p ** (e + math.floor(math.log(n) / math.log(p))) * term.evalf(subs={'x': p})))

    return result

def factor_combinations(factors):
    if not factors:
        yield ()
        return

    p, e = factors[0]
    for i in range(e + 1):
        for result in factor_combinations(factors[1:]):
            yield ((p, i),) + result

def product(factors):
    result = 1
    for p, e in factors:
        result *= p ** e

    return result

def divisors(n):
    factors = factor(n)
    for combination in factor_combinations(factors):
        #yield combination
        yield product(combination)


def lambdas():
    for k in range(2, 10):
        for n in range(2, 10):
            print "lambda(n=%d, k=%d) = %s" % (n, k, lambda_n(n, k))

def phis():
    for k in range(2, 10):
        for n in range(2, 10):
            print "phi(n=%d, k=%d) = %s" % (n, k, phi_n(n, k))

#phis()
phi = phi_n(15, 12)
print phi
print lambda_n(15, 12)
