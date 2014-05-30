#ifndef __MATRIX_H__
#define __MATRIX_H__
#include <cstdio>
#include <cstring>
#include <ctime>

#include <algorithm>
#include <vector>
#include <gmpxx.h>

#include "functions.h"

#define DEBUG

typedef struct {
    mpz_class lambda;
    unsigned int confidence;
    mpz_class phi;
} LambdaResult;

bool cmp_factors(const Factor &a, const Factor &b) {
    mpz_class t1, t2;
    mpz_pow_ui(t1.get_mpz_t(), a.prime.get_mpz_t(), a.exponent);
    mpz_pow_ui(t2.get_mpz_t(), b.prime.get_mpz_t(), b.exponent);
    return t1 > t2;
}

template <unsigned int N, unsigned int K>
class Matrix {
 public:
    inline Matrix() {
        m = new unsigned int[N * N];
    }

    inline ~Matrix() {
        if (m) {
            delete[] m;
        }
    }

    inline void operator=(const Matrix<N, K> &other) {
        memcpy(m, other.m, sizeof(unsigned int) * N * N);
    }

    inline void zero() {
        memset(m, 0, sizeof(unsigned int) * N * N);
    }

    inline void randomize() {
        for (unsigned int i = 0; i < N * N; ++i) {
            m[i] = ((unsigned int) rand() % (2 * 1000000)) / 1000000;
        }
    }

    inline bool next() {
        unsigned int i = 0;
        while (i < N * N) {
            m[i] = (m[i] + 1) % K;
            if (m[i] != 0) {
                return true;
            }
            ++i;
        }
        return false;
    }

    inline bool next_triangular() {
        unsigned int i = 0;
        while (i < N * N) {
            m[i] = (m[i] + 1) % K;
            if (m[i] != 0) {
                return true;
            }
            if (i % N == i / N) {
                // set diagonal element to 1 after it turned 0
                m[i] = 1;
            }
            ++i;
            while (i < N * N && i % N < i / N) {
                ++i;
            }
        }
        return false;
    }

    inline void make_identity() {
        zero();
        for (unsigned int i = 0; i < N * N; i += N + 1) {
            m[i] = 1;
        }
    }

    inline bool is_identity() const {
        if (N == 1) {
            return m[0] == 1;
        } else if (N == 2) {
            return m[0] == 1 && m[1] == 0 &&
                   m[2] == 0 && m[3] == 1;
        } else if (N == 3) {
            return m[0] == 1 && m[1] == 0 && m[2] == 0 &&
                   m[3] == 0 && m[4] == 1 && m[5] == 0 &&
                   m[6] == 0 && m[7] == 0 && m[8] == 1;
        } else if (N == 4) {
            return m[0] == 1 && m[1] == 0 && m[2] == 0 && m[3] == 0 &&
                   m[4] == 0 && m[5] == 1 && m[6] == 0 && m[7] == 0 &&
                   m[8] == 0 && m[9] == 0 && m[10] == 1 && m[11] == 0 &&
                   m[12] == 0 && m[13] == 0 && m[14] == 0 && m[15] == 1;
        } else if (N == 5) {
            return m[0] == 1 && m[1] == 0 && m[2] == 0 && m[3] == 0 && m[4] == 0 &&
                   m[5] == 0 && m[6] == 1 && m[7] == 0 && m[8] == 0 && m[9] == 0 &&
                   m[10] == 0 && m[11] == 0 && m[12] == 1 && m[13] == 0 && m[14] == 0 &&
                   m[15] == 0 && m[16] == 0 && m[17] == 0 && m[18] == 1 && m[19] == 0 &&
                   m[20] == 0 && m[21] == 0 && m[22] == 0 && m[23] == 0 && m[24] == 1;
        } else if (N == 6) {
            return m[0] == 1 && m[1] == 0 && m[2] == 0 && m[3] == 0 && m[4] == 0 && m[5] == 0 &&
                   m[6] == 0 && m[7] == 1 && m[8] == 0 && m[9] == 0 && m[10] == 0 && m[11] == 0 &&
                   m[12] == 0 && m[13] == 0 && m[14] == 1 && m[15] == 0 && m[16] == 0 && m[17] == 0 &&
                   m[18] == 0 && m[19] == 0 && m[20] == 0 && m[21] == 1 && m[22] == 0 && m[23] == 0 &&
                   m[24] == 0 && m[25] == 0 && m[26] == 0 && m[27] == 0 && m[28] == 1 && m[29] == 0 &&
                   m[30] == 0 && m[31] == 0 && m[32] == 0 && m[33] == 0 && m[34] == 0 && m[35] == 1;
        }
        unsigned int d = 0;
        for (unsigned int i = 0; i < N * N; ++i) {
            if (i == d) {
                if (m[i] != 1) {
                    return false;
                }
                d += N + 1;
            } else if (m[i] != 0) {
                return false;
            }
        }

        return true;
    }

    inline bool is_zero() const {
        for (unsigned int i = 0; i < N * N; ++i) {
            if (m[i] != 0) {
                return false;
            }
        }
        return true;
    }

    inline unsigned int &at(unsigned int row, unsigned int column) {
        return m[row * N + column];
    }

    inline const unsigned int &at(unsigned int row, unsigned int column) const {
        return m[row * N + column];
    }

    inline void mult(const Matrix<N, K> &other, Matrix<N, K> &result) const {
        for (unsigned int r = 0; r < N; ++r) {
            for (unsigned int c = 0; c < N; ++c) {
                unsigned int sum = 0;
                for (unsigned int i = 0; i < N; ++i) {
                    sum += at(r, i) * other.at(i, c);
                }
                result.at(r, c) = sum % K;
            }
        }
#ifdef DEBUG
        ++num_multiplications;
#endif
    }

    inline void mult(const Matrix<N, K> &other) {
        for (unsigned int r = 0; r < N; ++r) {
            for (unsigned int c = 0; c < N; ++c) {
                unsigned int sum = 0;
                for (unsigned int i = 0; i < N; ++i) {
                    sum += at(r, i) * other.at(i, c);
                }
                buffer[r * N + c] = sum % K;
            }
        }
        //memcpy(m, buffer, sizeof(unsigned int) * N * N);
        unsigned int *tmp = buffer;
        buffer = m;
        m = tmp;
#ifdef DEBUG
        ++num_multiplications;
#endif
    }

    inline void square() {
        for (unsigned int r = 0; r < N; ++r) {
            for (unsigned int c = 0; c < N; ++c) {
                unsigned int sum = 0;
                for (unsigned int i = 0; i < N; ++i) {
                    sum += at(r, i) * at(i, c);
                }
                buffer[r * N + c] = sum % K;
            }
        }
        //memcpy(m, buffer, sizeof(unsigned int) * N * N);
        unsigned int *tmp = buffer;
        buffer = m;
        m = tmp;
#ifdef DEBUG
        ++num_multiplications;
#endif
    }

    void power(const mpz_class &e, Matrix<N, K> &result) const {
        if (e == 0) {
            result.make_identity();
            return;
        } else if (e == 1) {
            result = *this;
            return;
        }

        power(e / 2, result);
        //result.mult(result);
        result.square();
        if (e % 2 == 1) {
            result.mult(*this);
        }
    }

    void print() const {
        for (unsigned int r = 0; r < N; ++r) {
            for (unsigned int c = 0; c < N; ++c) {
                printf("%2d", at(r, c));
                if (c == N - 1) {
                    printf("\n");
                } else {
                    printf(" ");
                }
            }
        }
    }

    static unsigned int search_for_prime_exponent(const Matrix<N, K> &base, const Factor &factor) {
        Matrix<N, K> tmp;
        unsigned int result = 0;
        mpz_class exp = 1;
        if (base.is_identity()) {
            return 0;
        }
        // this could be done with binary search, but in practice it seems like exp is usually low,
        // so naively searching from the start seems to be the fastest way
        for (unsigned int i = 0; i < factor.exponent; ++i) {
            exp *= factor.prime;
            ++result;
            base.power(exp, tmp);
            if (tmp.is_identity()) {
                return result;
            }
        }
        return result;
    }

    mpz_class get_minimal_exponent(const mpz_class &phi, const Factors &factors, Factors *left_factors) {
        mpz_class exponent = phi;
        Matrix<N, K> tmp;
        Factors new_factors = factors;

        for (unsigned int i = 0; i < factors.size(); ++i) {
            for (unsigned int j = 0; j < factors[i].exponent; ++j) {
                exponent /= factors[i].prime;
            }
            power(exponent, tmp);
            unsigned int exp = search_for_prime_exponent(tmp, factors[i]);
            //printf("%s : %d : %d\n", factors[i].prime.get_str().c_str(), factors[i].exponent, exp);
            new_factors[i].exponent -= exp;
            for (unsigned int j = 0; j < exp; ++j) {
                exponent *= factors[i].prime;
            }
        }

        if (left_factors) {
            *left_factors = new_factors;
        }
        return exponent;
    }

    static LambdaResult find_probable_lambda() {
        mpz_class phi = phi_n(N, K);
        Factors factors = factorize(phi);
        // go backwards, so large factors are removed early on
        sort(factors.begin(), factors.end(), cmp_factors);

        mpz_class num_matrices = 1;
        for (unsigned int i = 0; i < N * N; ++i) {
            num_matrices *= K;
        }
        mpz_class result = 1;
        bool check_all = false;
        if (phi < 50000000) {
            check_all = true;
        }
        Matrix<N, K> m;
        m.zero();
        unsigned int z = 0;
        unsigned int s = 0;
        mpz_class c, max_c = 200000, ctmp;
        if (check_all) {
            max_c = phi;
        }
        Matrix<N, K> tmp, tmp2;
        unsigned int last_timestamp = time(NULL);
        while (true) {
            if (check_all) {
                if (!m.next()) {
                    break;
                }
            } else {
                if (c >= max_c) {
                    break;
                }
                m.randomize();
            }
            if (time(NULL) >= last_timestamp + 5) {
                ctmp = 10000 * c / max_c;
                printf("%.f%% c:%s, z:%d, s:%d\n", (float) ctmp.get_ui() / 100, c.get_str().c_str(), z, s);
                last_timestamp = time(NULL);
            }

            //printf("%d %d\n", c, z);
            //m.print();
            mpz_class e = 1;
            m.power(phi, tmp);
            if (tmp.is_zero()) {
                ++z;
                continue;
            } else if (!tmp.is_identity()) {
                ++s;
                // singular matrix, ignore
                continue;
            }

            m.power(result, tmp);
            if (tmp.is_identity()) {
                ++c;
                continue;
            }

            /*for (unsigned int i = 0; i < factors.size(); ++i) {
                printf("%s^%d, ", factors[i].prime.get_str().c_str(), factors[i].exponent);
            }
            printf("\n");*/
            e = tmp.get_minimal_exponent(phi / result, factors, &factors);
            /*for (unsigned int i = 0; i < factors.size(); ++i) {
                printf("%s^%d, ", factors[i].prime.get_str().c_str(), factors[i].exponent);
            }
            printf("\n");*/
            ++c;
            result *= e;
            printf("%s\n", result.get_str().c_str());
            unsigned int prime_factors = 0;
            for (unsigned int i = 0; i < factors.size(); ++i) {
                prime_factors += factors[i].exponent;
            }
            printf("rest: %d\n", prime_factors);
        }

        unsigned int confidence;
        if (check_all) {
            confidence = 0;
        } else {
            confidence = c.get_ui();
        }
        return {result, confidence, num_matrices};
    }

#ifdef DEBUG
    static unsigned long int num_multiplications;
#endif

 private:
    unsigned int *m;
    static unsigned int *buffer;
};

template <unsigned int N, unsigned int K>
unsigned int *Matrix<N, K>::buffer = new unsigned int[N * N];

#ifdef DEBUG
template <unsigned int N, unsigned int K>
unsigned long int Matrix<N, K>::num_multiplications = 0;
#endif

#endif
