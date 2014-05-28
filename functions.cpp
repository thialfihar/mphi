#include "functions.h"

Factors factorize(mpz_class n) {
    Factors result;
    mpz_class p = 2;
    while (p * p <= n) {
        unsigned int e = 0;
        while (n % p == 0) {
            ++e;
            n /= p;
        }
        if (e > 0) {
            result.push_back(Factor(p, e));
        }
        mpz_nextprime(p.get_mpz_t(), p.get_mpz_t());
    }
    if (n != 1) {
        result.push_back(Factor(n, 1));
    }
    return result;
}

vector<mpz_class> get_divisors(mpz_class n) {
    Factors factors = factorize(n);
    vector<mpz_class> result;
    result.push_back(1);
    for (Factor f : factors) {
        auto tmp_results = result;
        for (unsigned int i = 1; i <= f.second; ++i) {
            mpz_class tmp;
            mpz_pow_ui(tmp.get_mpz_t(), f.first.get_mpz_t(), i);
            for (auto v : tmp_results) {
                result.push_back(tmp * v);
            }
        }
    }

    sort(result.begin(), result.end());
    return result;
}

mpz_class phi_n(unsigned int n, unsigned int k) {
    Factors factors = factorize(k);
    mpz_class result = 1;
    for (unsigned int i = 1; i <= n; ++i) {
        for (Factor f : factors) {
            mpz_class tmp1, tmp2;
            mpz_pow_ui(tmp1.get_mpz_t(), f.first.get_mpz_t(), f.second * n);
            mpz_pow_ui(tmp2.get_mpz_t(), f.first.get_mpz_t(), f.second * n - i);
            result *= tmp1 - tmp2;
        }
    }

    return result;
}
