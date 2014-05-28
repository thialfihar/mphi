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

    inline bool is_identity() {
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

    inline bool is_zero() {
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

    inline void mult(const Matrix<N, K> &other, Matrix<N, K> &result) {
        if (N == 1) {
            result.m[0] = (m[0] * other.m[0]) % K;
            return;
        } else if (N == 2) {
            result.m[0] = (m[0] * other.m[0] + m[1] * other.m[2]) % K;
            result.m[1] = (m[0] * other.m[1] + m[1] * other.m[3]) % K;
            result.m[2] = (m[2] * other.m[0] + m[3] * other.m[2]) % K;
            result.m[3] = (m[2] * other.m[1] + m[3] * other.m[3]) % K;
            return;
        } else if (N == 3) {
            result.m[0] = (m[0] * other.m[0] + m[1] * other.m[3] + m[2] * other.m[6]) % K;
            result.m[1] = (m[0] * other.m[1] + m[1] * other.m[4] + m[2] * other.m[7]) % K;
            result.m[2] = (m[0] * other.m[2] + m[1] * other.m[5] + m[2] * other.m[8]) % K;
            result.m[3] = (m[3] * other.m[0] + m[4] * other.m[3] + m[5] * other.m[6]) % K;
            result.m[4] = (m[3] * other.m[1] + m[4] * other.m[4] + m[5] * other.m[7]) % K;
            result.m[5] = (m[3] * other.m[2] + m[4] * other.m[5] + m[5] * other.m[8]) % K;
            result.m[6] = (m[6] * other.m[0] + m[7] * other.m[3] + m[8] * other.m[6]) % K;
            result.m[7] = (m[6] * other.m[1] + m[7] * other.m[4] + m[8] * other.m[7]) % K;
            result.m[8] = (m[6] * other.m[2] + m[7] * other.m[5] + m[8] * other.m[8]) % K;
            return;
        } else if (N == 4) {
            result.m[0] = (m[0] * other.m[0] + m[1] * other.m[4] + m[2] * other.m[8] + m[3] * other.m[12]) % K;
            result.m[1] = (m[0] * other.m[1] + m[1] * other.m[5] + m[2] * other.m[9] + m[3] * other.m[13]) % K;
            result.m[2] = (m[0] * other.m[2] + m[1] * other.m[6] + m[2] * other.m[10] + m[3] * other.m[14]) % K;
            result.m[3] = (m[0] * other.m[3] + m[1] * other.m[7] + m[2] * other.m[11] + m[3] * other.m[15]) % K;
            result.m[4] = (m[4] * other.m[0] + m[5] * other.m[4] + m[6] * other.m[8] + m[7] * other.m[12]) % K;
            result.m[5] = (m[4] * other.m[1] + m[5] * other.m[5] + m[6] * other.m[9] + m[7] * other.m[13]) % K;
            result.m[6] = (m[4] * other.m[2] + m[5] * other.m[6] + m[6] * other.m[10] + m[7] * other.m[14]) % K;
            result.m[7] = (m[4] * other.m[3] + m[5] * other.m[7] + m[6] * other.m[11] + m[7] * other.m[15]) % K;
            result.m[8] = (m[8] * other.m[0] + m[9] * other.m[4] + m[10] * other.m[8] + m[11] * other.m[12]) % K;
            result.m[9] = (m[8] * other.m[1] + m[9] * other.m[5] + m[10] * other.m[9] + m[11] * other.m[13]) % K;
            result.m[10] = (m[8] * other.m[2] + m[9] * other.m[6] + m[10] * other.m[10] + m[11] * other.m[14]) % K;
            result.m[11] = (m[8] * other.m[3] + m[9] * other.m[7] + m[10] * other.m[11] + m[11] * other.m[15]) % K;
            result.m[12] = (m[12] * other.m[0] + m[13] * other.m[4] + m[14] * other.m[8] + m[15] * other.m[12]) % K;
            result.m[13] = (m[12] * other.m[1] + m[13] * other.m[5] + m[14] * other.m[9] + m[15] * other.m[13]) % K;
            result.m[14] = (m[12] * other.m[2] + m[13] * other.m[6] + m[14] * other.m[10] + m[15] * other.m[14]) % K;
            result.m[15] = (m[12] * other.m[3] + m[13] * other.m[7] + m[14] * other.m[11] + m[15] * other.m[15]) % K;
            return;
        }

        for (unsigned int r = 0; r < N; ++r) {
            for (unsigned int c = 0; c < N; ++c) {
                unsigned int sum = 0;
                for (unsigned int i = 0; i < N; ++i) {
                    sum += m[r * N + i] * other.m[i * N + c];
                }
                result.m[r * N + c] = sum % K;
            }
        }
    }

    inline void mult(const Matrix<N, K> &other) {
        if (N == 1) {
            buffer[0] = (m[0] * other.m[0]) % K;
            memcpy(m, buffer, sizeof(unsigned int) * N * N);
            return;
        } else if (N == 2) {
            buffer[0] = (m[0] * other.m[0] + m[1] * other.m[2]) % K;
            buffer[1] = (m[0] * other.m[1] + m[1] * other.m[3]) % K;
            buffer[2] = (m[2] * other.m[0] + m[3] * other.m[2]) % K;
            buffer[3] = (m[2] * other.m[1] + m[3] * other.m[3]) % K;
            memcpy(m, buffer, sizeof(unsigned int) * N * N);
            return;
        } else if (N == 3) {
            buffer[0] = (m[0] * other.m[0] + m[1] * other.m[3] + m[2] * other.m[6]) % K;
            buffer[1] = (m[0] * other.m[1] + m[1] * other.m[4] + m[2] * other.m[7]) % K;
            buffer[2] = (m[0] * other.m[2] + m[1] * other.m[5] + m[2] * other.m[8]) % K;
            buffer[3] = (m[3] * other.m[0] + m[4] * other.m[3] + m[5] * other.m[6]) % K;
            buffer[4] = (m[3] * other.m[1] + m[4] * other.m[4] + m[5] * other.m[7]) % K;
            buffer[5] = (m[3] * other.m[2] + m[4] * other.m[5] + m[5] * other.m[8]) % K;
            buffer[6] = (m[6] * other.m[0] + m[7] * other.m[3] + m[8] * other.m[6]) % K;
            buffer[7] = (m[6] * other.m[1] + m[7] * other.m[4] + m[8] * other.m[7]) % K;
            buffer[8] = (m[6] * other.m[2] + m[7] * other.m[5] + m[8] * other.m[8]) % K;
            memcpy(m, buffer, sizeof(unsigned int) * N * N);
            return;
        } else if (N == 4) {
            buffer[0] = (m[0] * other.m[0] + m[1] * other.m[4] + m[2] * other.m[8] + m[3] * other.m[12]) % K;
            buffer[1] = (m[0] * other.m[1] + m[1] * other.m[5] + m[2] * other.m[9] + m[3] * other.m[13]) % K;
            buffer[2] = (m[0] * other.m[2] + m[1] * other.m[6] + m[2] * other.m[10] + m[3] * other.m[14]) % K;
            buffer[3] = (m[0] * other.m[3] + m[1] * other.m[7] + m[2] * other.m[11] + m[3] * other.m[15]) % K;
            buffer[4] = (m[4] * other.m[0] + m[5] * other.m[4] + m[6] * other.m[8] + m[7] * other.m[12]) % K;
            buffer[5] = (m[4] * other.m[1] + m[5] * other.m[5] + m[6] * other.m[9] + m[7] * other.m[13]) % K;
            buffer[6] = (m[4] * other.m[2] + m[5] * other.m[6] + m[6] * other.m[10] + m[7] * other.m[14]) % K;
            buffer[7] = (m[4] * other.m[3] + m[5] * other.m[7] + m[6] * other.m[11] + m[7] * other.m[15]) % K;
            buffer[8] = (m[8] * other.m[0] + m[9] * other.m[4] + m[10] * other.m[8] + m[11] * other.m[12]) % K;
            buffer[9] = (m[8] * other.m[1] + m[9] * other.m[5] + m[10] * other.m[9] + m[11] * other.m[13]) % K;
            buffer[10] = (m[8] * other.m[2] + m[9] * other.m[6] + m[10] * other.m[10] + m[11] * other.m[14]) % K;
            buffer[11] = (m[8] * other.m[3] + m[9] * other.m[7] + m[10] * other.m[11] + m[11] * other.m[15]) % K;
            buffer[12] = (m[12] * other.m[0] + m[13] * other.m[4] + m[14] * other.m[8] + m[15] * other.m[12]) % K;
            buffer[13] = (m[12] * other.m[1] + m[13] * other.m[5] + m[14] * other.m[9] + m[15] * other.m[13]) % K;
            buffer[14] = (m[12] * other.m[2] + m[13] * other.m[6] + m[14] * other.m[10] + m[15] * other.m[14]) % K;
            buffer[15] = (m[12] * other.m[3] + m[13] * other.m[7] + m[14] * other.m[11] + m[15] * other.m[15]) % K;
            memcpy(m, buffer, sizeof(unsigned int) * N * N);
            return;
        }

        for (unsigned int r = 0; r < N; ++r) {
            for (unsigned int c = 0; c < N; ++c) {
                unsigned int sum = 0;
                for (unsigned int i = 0; i < N; ++i) {
                    sum += m[r * N + i] * other.m[i * N + c];
                }
                buffer[r * N + c] = sum % K;
            }
        }
        memcpy(m, buffer, sizeof(unsigned int) * N * N);
    }

    void power(const mpz_class &e, Matrix<N, K> &result) {
        if (e == 0) {
            result.make_identity();
            return;
        } else if (e == 1) {
            result = *this;
            return;
        }

        power(e / 2, result);
        result.mult(result);
        if (e % 2 == 1) {
            result.mult(*this);
        }
    }

    void print() {
        for (unsigned int r = 0; r < N; ++r) {
            for (unsigned int c = 0; c < N; ++c) {
                printf("%8d", at(r, c));
                if (c == N - 1) {
                    printf("\n");
                } else {
                    printf("\t");
                }
            }
        }
    }

    //static mpz_class find_lambda(unsigned int K);
    //static mpz_class find_probable_lambda(unsigned int K);

    static mpz_class find_lambda() {
        mpz_class phi = phi_n(N, K);
        vector<mpz_class> candidates = get_divisors(phi);
        mpz_class result = 1;
        if (phi > 100000000) {
            //return 0;
        }

        Matrix<N, K> m;
        m.zero();
        // for triangular approach
        //m.make_identity();
        unsigned int c = 0;
        while (m.next() && c < 20000) {
            //m.print();
            //++c;
            Matrix<N, K> tmp;
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

    static mpz_class find_probable_lambda() {
        mpz_class phi = phi_n(N, K);
        vector<mpz_class> candidates = get_divisors(phi);
        mpz_class result = 1;
        if (phi < 20000) {
            return find_lambda();
        }

        Matrix<N, K> m;
        unsigned int c = 0;
        unsigned int z = 0;
        unsigned int max_c = 100;
        while (c < max_c) {
            m.randomize();
            //printf("%d %d\n", c, z);
            //m.print();
            Matrix<N, K> tmp;
            mpz_class e = 1;
            m.power(result, tmp);
            if (tmp.is_zero()) {
                ++z;
                if (z % 10 == 0) {
                    printf("z: %d\n", z);
                }
                continue;
            } else if (tmp.is_identity()) {
                ++c;
                if (c % 10 == 0) {
                    printf("%.f%%\n", 100.0 * c / max_c);
                }
                continue;
            }
            for (unsigned int i = 0; i < candidates.size(); ++i) {
                //printf("%d/%d\n", i, candidates.size());
                m.power(candidates[i], tmp);
                if (tmp.is_zero()) {
                    ++z;
                    if (z % 10 == 0) {
                        printf("z: %d\n", z);
                    }
                    e = 1;
                    //printf("break at %d\n", i);
                    break;
                } else if (tmp.is_identity()) {
                    e = candidates[i];
                    ++c;
                    if (c % 10 == 0) {
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

 private:
    unsigned int m[N * N];
    static unsigned int buffer[N * N];
};

template <unsigned int N, unsigned int K>
unsigned int Matrix<N, K>::buffer[N * N];

#endif
