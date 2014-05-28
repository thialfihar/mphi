#ifndef __MATRIX_H__
#define __MATRIX_H__
#include <cstdio>
#include <cstring>
#include <ctime>

#include <algorithm>
#include <vector>
#include <gmpxx.h>

#include "functions.h"

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
            m[i] = ((unsigned int) rand() % (K * 1000000)) / 1000000;
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

    mpz_class get_minimal_exponent(const mpz_class &phi, const Factors &factors) {
        mpz_class exponent = phi;
        Matrix<N, K> tmp;
        for (unsigned int i = 0; i < factors.size(); ++i) {
            mpz_class e = exponent;
            for (unsigned int j = 0; j < factors[i].second; ++j) {
                e /= factors[i].first;
                power(e, tmp);
                if (!tmp.is_identity()) {
                    break;
                }
                exponent = e;
            }
        }

        return exponent;
    }

    static mpz_class find_lambda() {
        mpz_class phi = phi_n(N, K);
        Factors factors = factorize(phi);
        vector<mpz_class> candidates = get_divisors(phi);
        mpz_class result = 1;

        Matrix<N, K> m;
        m.zero();
        // for triangular approach
        //m.make_identity();
        unsigned int c = 0;
        Matrix<N, K> tmp;
        while (m.next()) {
            //m.print();
            //++c;
            mpz_class e = 1;
            m.power(phi, tmp);
            if (tmp.is_zero()) {
                continue;
            } else if (!tmp.is_identity()) {
                // singular matrix, ignore
                continue;
            }

            m.power(result, tmp);
            if (tmp.is_zero()) {
                continue;
            } else if (tmp.is_identity()) {
                continue;
            }
            e = m.get_minimal_exponent(phi, factors);

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

    static mpz_class find_probable_lambda() {
        mpz_class phi = phi_n(N, K);
        Factors factors = factorize(phi);
        vector<mpz_class> candidates = get_divisors(phi);
        //printf("phi = %s, first = %s, last = %s\n", phi.get_str().c_str(), candidates[0].get_str().c_str(), candidates[candidates.size() - 1].get_str().c_str());
        candidates.erase(candidates.begin());
        candidates.pop_back();
        mpz_class result = 1;
        if (phi < 2000000) {
            return find_lambda();
        }

        Matrix<N, K> m;
        unsigned int c = 0;
        unsigned int z = 0;
        unsigned int s = 0;
        unsigned int max_c = 1000000;
        Matrix<N, K> tmp;
        while (c < max_c) {
            m.randomize();
            //printf("%d %d\n", c, z);
            //m.print();
            mpz_class e = 1;
            m.power(phi, tmp);
            if (tmp.is_zero()) {
                ++z;
                if (z % 100 == 0) {
                    printf("z: %d\n", z);
                }
                continue;
            } else if (!tmp.is_identity()) {
                ++s;
                // singular matrix, ignore
                continue;
            }

            m.power(result, tmp);
            if (tmp.is_identity()) {
                ++c;
                if (c % 100000 == 0) {
                    printf("%.f%%\n", 100.0 * c / max_c);
                }
                continue;
            }

            e = m.get_minimal_exponent(phi, factors);

            if (e > 1) {
                ++c;
                if (c % 100000 == 0) {
                    printf("%.f%%\n", 100.0 * c / max_c);
                }
                mpz_lcm(result.get_mpz_t(), result.get_mpz_t(), e.get_mpz_t());
                /*auto it = candidates.begin();
                while (it != candidates.end()) {
                    if (result % *it == 0) {
                        candidates.erase(it);
                    } else {
                        ++it;
                    }
                }
                printf("%d c: %d - %s\n", c, (unsigned int)candidates.size(), result.get_str().c_str());*/
            }
        }

        return result;
    }

 private:
    unsigned int *m;
    static unsigned int *buffer;
};

template <unsigned int N, unsigned int K>
unsigned int *Matrix<N, K>::buffer = new unsigned int[N * N];

#endif
