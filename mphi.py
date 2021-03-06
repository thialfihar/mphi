#!/usr/bin/env python3
import json
import locale
import math
import subprocess
from sympy import symbols, lcm

encoding = locale.getdefaultlocale()[1]

class Matrix:
    def __init__(self, n, k, id=None):
        self.n = n
        self.k = k
        self.m = [0 for i in range(n * n)]
        if id:
            i = 0
            while id:
                self.m[i] = id % k
                id //= k
                i += 1

    def get_id(self):
        id = 0
        for i in range(self.n * self.n - 1, -1, -1):
            id *= self.k
            id += self.m[i]

        return id

    def scalar_mult(self, s):
        s = ((s % self.k) + self.k) % self.k
        for i in range(self.n * self.n):
            self.m[i] = (self.m[i] * s) % self.k

    def print(self):
        for i in range(self.n):
            for j in range(self.n):
                print('%d ' % self.m[i * self.n + j], end='')
            print()

def factor(n):
    factors = {}
    for p in subprocess.check_output(["factor", str(n)]).decode(encoding).split(':')[1].strip().split(' '):
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
            print("lambda(n=%d, k=%d) = %s" % (n, k, lambda_n(n, k)))

def phis():
    for k in range(2, 10):
        for n in range(2, 10):
            print("phi(n=%d, k=%d) = %s" % (n, k, phi_n(n, k)))

def analyse_map(n, k):
    data = json.load(open('map_%d_%d.json' % (n, k), 'r'))
    factors = factor(k)
    phi = phi_n(n, k)
    x = symbols('p')
    terms = []
    for i in range(1, n + 1):
        terms = [x ** i - 1] + terms

    exponents = map(str, sorted(map(int, data['map'].keys())))
    for exponent in exponents:
        num_matrices = data['map'][exponent]
        exponent = int(exponent)
        print(exponent)
        tmp = exponent
        tmp2 = 0
        while tmp > 1 and tmp % p == 0:
            tmp /= p
            tmp2 += 1
        if tmp2 > 0:
            print("    = %s * %s [%s]" % (tmp, exponent / tmp, 'p**%d' % tmp2))

        for p, e in factors:
            for term in terms:
                f = int(term.evalf(subs={'p': p}))
                if exponent % f == 0:
                    print("    = %s * %s [%s]" % (exponent / f, f, term))

    counts = {}
    squares = set()
    for cycle in data['cycles']:
        if len(cycle) == 2:
            squares.add(cycle[0])
        if len(cycle) % 2 == 0:
            squares.add(cycle[len(cycle) // 2 - 1])
        for mid in cycle:
            counts[mid] = counts.get(mid, 0) + 1

    keys = list(counts.keys())
    keys.sort(key=lambda x: counts[x])
    for mid in keys:
        if mid not in squares:
            continue
        print(mid, counts[mid], mid in squares)
        Matrix(n, k, mid).print()
        print()



analyse_map(3, 3)

m = Matrix(3, 3, 14662)
#m.m = [2, 0, 0, 0, 2, 0, 0, 0, 2]
m.print()
print(m.get_id())
m.scalar_mult(-1)
m.print()
print(m.get_id())
