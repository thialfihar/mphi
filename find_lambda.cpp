#include <cstdio>
#include <cstring>

#include <algorithm>
#include <vector>
#include <gmpxx.h>

#include "Matrix.h"

using std::pair;
using std::vector;

typedef pair<mpz_class, unsigned int> Factor;
typedef vector<Factor> Factors;

unsigned int Matrix::buffer[MAX_MATRIX_SIZE * MAX_MATRIX_SIZE];

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

mpz_class find_lambda(unsigned int n, unsigned int k) {
    mpz_class phi = phi_n(n, k);
    vector<mpz_class> candidates = get_divisors(phi);
    mpz_class result = 1;
    if (phi > 100000000) {
        //return 0;
    }

    Matrix m(n, k);
    m.zero();
    // for triangular approach
    //m.make_identity();
    unsigned int c = 0;
    while (m.next() && c < 20000) {
        //m.print();
        //++c;
        Matrix tmp(n, k);
        mpz_class e = 1;
        m.power(result, tmp);
        if (tmp.is_zero()) {
            continue;
        } else if (tmp.is_identity()) {
            continue;
        }
        for (unsigned int i = 0; i < candidates.size(); ++i) {
            m.power(candidates[i], tmp);
            if (tmp.is_zero()) {
                e = 1;
                //printf("break at %d\n", i);
                break;
            } else if (tmp.is_identity()) {
                e = candidates[i];
                break;
            }
        }
        if (e > 1) {
            mpz_lcm(result.get_mpz_t(), result.get_mpz_t(), e.get_mpz_t());
            auto it = candidates.begin();
            while (it != candidates.end()) {
                if (result % *it == 0) {
                    candidates.erase(it);
                } else {
                    ++it;
                }
            }
            //printf("c: %d - %s\n", (unsigned int)candidates.size(), result.get_str().c_str());
        }
    }

    return result;
}

mpz_class find_probable_lambda(unsigned int n, unsigned int k) {
    mpz_class phi = phi_n(n, k);
    vector<mpz_class> candidates = get_divisors(phi);
    mpz_class result = 1;
    if (phi < 20000) {
        return find_lambda(n, k);
    }

    Matrix m(n, k);
    unsigned int c = 0;
    unsigned int z = 0;
    unsigned int max_c = 500;
    while (c < max_c) {
        m.randomize();
        //printf("%d %d\n", c, z);
        //m.print();
        Matrix tmp(n, k);
        mpz_class e = 1;
        m.power(result, tmp);
        if (tmp.is_zero()) {
            ++z;
            if (z % 1000 == 0) {
                printf("z: %d\n", z);
            }
            continue;
        } else if (tmp.is_identity()) {
            ++c;
            if (c % 1000 == 0) {
                printf("%.f%%\n", 100.0 * c / max_c);
            }
            continue;
        }
        for (unsigned int i = 0; i < candidates.size(); ++i) {
            //printf("%d/%d\n", i, candidates.size());
            m.power(candidates[i], tmp);
            if (tmp.is_zero()) {
                ++z;
                if (z % 1000 == 0) {
                    printf("z: %d\n", z);
                }
                e = 1;
                //printf("break at %d\n", i);
                break;
            } else if (tmp.is_identity()) {
                e = candidates[i];
                ++c;
                if (c % 1000 == 0) {
                    printf("%.f%%\n", 100.0 * c / max_c);
                }
                break;
            }
        }
        if (e > 1) {
            mpz_lcm(result.get_mpz_t(), result.get_mpz_t(), e.get_mpz_t());
            auto it = candidates.begin();
            while (it != candidates.end()) {
                if (result % *it == 0) {
                    candidates.erase(it);
                } else {
                    ++it;
                }
            }
            printf("%d c: %d - %s\n", c, (unsigned int)candidates.size(), result.get_str().c_str());
        }
    }

    return result;
}

void test_stuff() {
    Matrix m1(3, 7);
    Matrix m2(3, 7);
    Matrix result(3, 7);

    m1.make_identity();
    m2.make_identity();

    m1.mult(m2, result);
    m2.power(134, result);
    result.print();
    if (result.is_identity()) {
        printf("yup\n");
        result.print();
    }

    Factors factors = factorize(193634213232);
    for (Factor f : factors) {
        printf("%s^%d\n", f.first.get_str().c_str(), f.second);
    }

    printf("%s\n", phi_n(7, 2).get_str().c_str());

    printf("%d\n", (unsigned int)get_divisors(163849992929280).size());
    for (auto v : get_divisors(163849992929280)) {
        printf("%s\n", v.get_str().c_str());
    }
}

int main(int argc, char *argv[]) {
    srand(time(nullptr));
    //test_stuff(); return 0;

    for (unsigned int k = 6; k <= 6; ++k) {
        for (unsigned int n = 5; n <= 5; ++n) {
            printf("lambda(n=%d, k=%d) = %s\n", n, k, find_probable_lambda(n, k).get_str().c_str());
            fflush(stdout);
        }
    }
    return 0;
}
