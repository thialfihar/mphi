#ifndef __GENERAL_MATRIX_H__
#define __GENERAL_MATRIX_H__
#include <cstdio>
#include <cstring>
#include <ctime>

#include <algorithm>
#include <vector>
#include <gmpxx.h>

#include "functions.h"

namespace gm {

class Matrix {
 public:
    inline Matrix(unsigned int n, unsigned int k)
        :n(n), k(k) {
    }

    inline void operator=(const Matrix &other) {
        memcpy(m, other.m, sizeof(unsigned int) * n * n);
    }

    inline void zero() {
        memset(m, 0, sizeof(unsigned int) * n * n);
    }

    inline void randomize() {
        for (unsigned int i = 0; i < n * n; ++i) {
            m[i] = ((unsigned int) rand() % (k * 1000000)) / 1000000;
        }
    }

    inline bool next() {
        unsigned int i = 0;
        while (i < n * n) {
            m[i] = (m[i] + 1) % k;
            if (m[i] != 0) {
                return true;
            }
            ++i;
        }
        return false;
    }

    inline bool next_triangular() {
        unsigned int i = 0;
        while (i < n * n) {
            m[i] = (m[i] + 1) % k;
            if (m[i] != 0) {
                return true;
            }
            if (i % n == i / n) {
                // set diagonal element to 1 after it turned 0
                m[i] = 1;
            }
            ++i;
            while (i < n * n && i % n < i / n) {
                ++i;
            }
        }
        return false;
    }

    inline void make_identity() {
        zero();
        for (unsigned int i = 0; i < n * n; i += n + 1) {
            m[i] = 1;
        }
    }

    inline bool is_identity() {
        if (n == 1) {
            return m[0] == 1;
        } else if (n == 2) {
            return m[0] == 1 && m[1] == 0 &&
                   m[2] == 0 && m[3] == 1;
        } else if (n == 3) {
            return m[0] == 1 && m[1] == 0 && m[2] == 0 &&
                   m[3] == 0 && m[4] == 1 && m[5] == 0 &&
                   m[6] == 0 && m[7] == 0 && m[8] == 1;
        } else if (n == 4) {
            return m[0] == 1 && m[1] == 0 && m[2] == 0 && m[3] == 0 &&
                   m[4] == 0 && m[5] == 1 && m[6] == 0 && m[7] == 0 &&
                   m[8] == 0 && m[9] == 0 && m[10] == 1 && m[11] == 0 &&
                   m[12] == 0 && m[13] == 0 && m[14] == 0 && m[15] == 1;
        } else if (n == 5) {
            return m[0] == 1 && m[1] == 0 && m[2] == 0 && m[3] == 0 && m[4] == 0 &&
                   m[5] == 0 && m[6] == 1 && m[7] == 0 && m[8] == 0 && m[9] == 0 &&
                   m[10] == 0 && m[11] == 0 && m[12] == 1 && m[13] == 0 && m[14] == 0 &&
                   m[15] == 0 && m[16] == 0 && m[17] == 0 && m[18] == 1 && m[19] == 0 &&
                   m[20] == 0 && m[21] == 0 && m[22] == 0 && m[23] == 0 && m[24] == 1;
        } else if (n == 6) {
            return m[0] == 1 && m[1] == 0 && m[2] == 0 && m[3] == 0 && m[4] == 0 && m[5] == 0 &&
                   m[6] == 0 && m[7] == 1 && m[8] == 0 && m[9] == 0 && m[10] == 0 && m[11] == 0 &&
                   m[12] == 0 && m[13] == 0 && m[14] == 1 && m[15] == 0 && m[16] == 0 && m[17] == 0 &&
                   m[18] == 0 && m[19] == 0 && m[20] == 0 && m[21] == 1 && m[22] == 0 && m[23] == 0 &&
                   m[24] == 0 && m[25] == 0 && m[26] == 0 && m[27] == 0 && m[28] == 1 && m[29] == 0 &&
                   m[30] == 0 && m[31] == 0 && m[32] == 0 && m[33] == 0 && m[34] == 0 && m[35] == 1;
        }
        unsigned int d = 0;
        for (unsigned int i = 0; i < n * n; ++i) {
            if (i == d) {
                if (m[i] != 1) {
                    return false;
                }
                d += n + 1;
            } else if (m[i] != 0) {
                return false;
            }
        }

        return true;
    }

    inline bool is_zero() {
        for (unsigned int i = 0; i < n * n; ++i) {
            if (m[i] != 0) {
                return false;
            }
        }
        return true;
    }

    inline unsigned int &at(unsigned int row, unsigned int column) {
        return m[row * n + column];
    }

    inline void mult(Matrix &other, Matrix &result) {
        result.zero();

        if (n == 1) {
            result.m[0] = (m[0] * other.m[0]) % k;
            return;
        } else if (n == 2) {
            result.m[0] = (m[0] * other.m[0] + m[1] * other.m[2]) % k;
            result.m[1] = (m[0] * other.m[1] + m[1] * other.m[3]) % k;
            result.m[2] = (m[2] * other.m[0] + m[3] * other.m[2]) % k;
            result.m[3] = (m[2] * other.m[1] + m[3] * other.m[3]) % k;
            return;
        } else if (n == 3) {
            result.m[0] = (m[0] * other.m[0] + m[1] * other.m[3] + m[2] * other.m[6]) % k;
            result.m[1] = (m[0] * other.m[1] + m[1] * other.m[4] + m[2] * other.m[7]) % k;
            result.m[2] = (m[0] * other.m[2] + m[1] * other.m[5] + m[2] * other.m[8]) % k;
            result.m[3] = (m[3] * other.m[0] + m[4] * other.m[3] + m[5] * other.m[6]) % k;
            result.m[4] = (m[3] * other.m[1] + m[4] * other.m[4] + m[5] * other.m[7]) % k;
            result.m[5] = (m[3] * other.m[2] + m[4] * other.m[5] + m[5] * other.m[8]) % k;
            result.m[6] = (m[6] * other.m[0] + m[7] * other.m[3] + m[8] * other.m[6]) % k;
            result.m[7] = (m[6] * other.m[1] + m[7] * other.m[4] + m[8] * other.m[7]) % k;
            result.m[8] = (m[6] * other.m[2] + m[7] * other.m[5] + m[8] * other.m[8]) % k;
            return;
        } else if (n == 4) {
            result.m[0] = (m[0] * other.m[0] + m[1] * other.m[4] + m[2] * other.m[8] + m[3] * other.m[12]) % k;
            result.m[1] = (m[0] * other.m[1] + m[1] * other.m[5] + m[2] * other.m[9] + m[3] * other.m[13]) % k;
            result.m[2] = (m[0] * other.m[2] + m[1] * other.m[6] + m[2] * other.m[10] + m[3] * other.m[14]) % k;
            result.m[3] = (m[0] * other.m[3] + m[1] * other.m[7] + m[2] * other.m[11] + m[3] * other.m[15]) % k;
            result.m[4] = (m[4] * other.m[0] + m[5] * other.m[4] + m[6] * other.m[8] + m[7] * other.m[12]) % k;
            result.m[5] = (m[4] * other.m[1] + m[5] * other.m[5] + m[6] * other.m[9] + m[7] * other.m[13]) % k;
            result.m[6] = (m[4] * other.m[2] + m[5] * other.m[6] + m[6] * other.m[10] + m[7] * other.m[14]) % k;
            result.m[7] = (m[4] * other.m[3] + m[5] * other.m[7] + m[6] * other.m[11] + m[7] * other.m[15]) % k;
            result.m[8] = (m[8] * other.m[0] + m[9] * other.m[4] + m[10] * other.m[8] + m[11] * other.m[12]) % k;
            result.m[9] = (m[8] * other.m[1] + m[9] * other.m[5] + m[10] * other.m[9] + m[11] * other.m[13]) % k;
            result.m[10] = (m[8] * other.m[2] + m[9] * other.m[6] + m[10] * other.m[10] + m[11] * other.m[14]) % k;
            result.m[11] = (m[8] * other.m[3] + m[9] * other.m[7] + m[10] * other.m[11] + m[11] * other.m[15]) % k;
            result.m[12] = (m[12] * other.m[0] + m[13] * other.m[4] + m[14] * other.m[8] + m[15] * other.m[12]) % k;
            result.m[13] = (m[12] * other.m[1] + m[13] * other.m[5] + m[14] * other.m[9] + m[15] * other.m[13]) % k;
            result.m[14] = (m[12] * other.m[2] + m[13] * other.m[6] + m[14] * other.m[10] + m[15] * other.m[14]) % k;
            result.m[15] = (m[12] * other.m[3] + m[13] * other.m[7] + m[14] * other.m[11] + m[15] * other.m[15]) % k;
            return;
        }

        for (unsigned int r = 0; r < n; ++r) {
            for (unsigned int c = 0; c < n; ++c) {
                unsigned int sum = 0;
                for (unsigned int i = 0; i < n; ++i) {
                    sum += m[r * n + i] * other.m[i * n + c];
                }
                result.m[r * 10 + c] = sum % k;
            }
        }
    }

    inline void mult(Matrix &other) {
        if (n == 1) {
            buffer[0] = (m[0] * other.m[0]) % k;
            memcpy(m, buffer, sizeof(unsigned int) * n * n);
            return;
        } else if (n == 2) {
            buffer[0] = (m[0] * other.m[0] + m[1] * other.m[2]) % k;
            buffer[1] = (m[0] * other.m[1] + m[1] * other.m[3]) % k;
            buffer[2] = (m[2] * other.m[0] + m[3] * other.m[2]) % k;
            buffer[3] = (m[2] * other.m[1] + m[3] * other.m[3]) % k;
            memcpy(m, buffer, sizeof(unsigned int) * n * n);
            return;
        } else if (n == 3) {
            buffer[0] = (m[0] * other.m[0] + m[1] * other.m[3] + m[2] * other.m[6]) % k;
            buffer[1] = (m[0] * other.m[1] + m[1] * other.m[4] + m[2] * other.m[7]) % k;
            buffer[2] = (m[0] * other.m[2] + m[1] * other.m[5] + m[2] * other.m[8]) % k;
            buffer[3] = (m[3] * other.m[0] + m[4] * other.m[3] + m[5] * other.m[6]) % k;
            buffer[4] = (m[3] * other.m[1] + m[4] * other.m[4] + m[5] * other.m[7]) % k;
            buffer[5] = (m[3] * other.m[2] + m[4] * other.m[5] + m[5] * other.m[8]) % k;
            buffer[6] = (m[6] * other.m[0] + m[7] * other.m[3] + m[8] * other.m[6]) % k;
            buffer[7] = (m[6] * other.m[1] + m[7] * other.m[4] + m[8] * other.m[7]) % k;
            buffer[8] = (m[6] * other.m[2] + m[7] * other.m[5] + m[8] * other.m[8]) % k;
            memcpy(m, buffer, sizeof(unsigned int) * n * n);
            return;
        } else if (n == 4) {
            buffer[0] = (m[0] * other.m[0] + m[1] * other.m[4] + m[2] * other.m[8] + m[3] * other.m[12]) % k;
            buffer[1] = (m[0] * other.m[1] + m[1] * other.m[5] + m[2] * other.m[9] + m[3] * other.m[13]) % k;
            buffer[2] = (m[0] * other.m[2] + m[1] * other.m[6] + m[2] * other.m[10] + m[3] * other.m[14]) % k;
            buffer[3] = (m[0] * other.m[3] + m[1] * other.m[7] + m[2] * other.m[11] + m[3] * other.m[15]) % k;
            buffer[4] = (m[4] * other.m[0] + m[5] * other.m[4] + m[6] * other.m[8] + m[7] * other.m[12]) % k;
            buffer[5] = (m[4] * other.m[1] + m[5] * other.m[5] + m[6] * other.m[9] + m[7] * other.m[13]) % k;
            buffer[6] = (m[4] * other.m[2] + m[5] * other.m[6] + m[6] * other.m[10] + m[7] * other.m[14]) % k;
            buffer[7] = (m[4] * other.m[3] + m[5] * other.m[7] + m[6] * other.m[11] + m[7] * other.m[15]) % k;
            buffer[8] = (m[8] * other.m[0] + m[9] * other.m[4] + m[10] * other.m[8] + m[11] * other.m[12]) % k;
            buffer[9] = (m[8] * other.m[1] + m[9] * other.m[5] + m[10] * other.m[9] + m[11] * other.m[13]) % k;
            buffer[10] = (m[8] * other.m[2] + m[9] * other.m[6] + m[10] * other.m[10] + m[11] * other.m[14]) % k;
            buffer[11] = (m[8] * other.m[3] + m[9] * other.m[7] + m[10] * other.m[11] + m[11] * other.m[15]) % k;
            buffer[12] = (m[12] * other.m[0] + m[13] * other.m[4] + m[14] * other.m[8] + m[15] * other.m[12]) % k;
            buffer[13] = (m[12] * other.m[1] + m[13] * other.m[5] + m[14] * other.m[9] + m[15] * other.m[13]) % k;
            buffer[14] = (m[12] * other.m[2] + m[13] * other.m[6] + m[14] * other.m[10] + m[15] * other.m[14]) % k;
            buffer[15] = (m[12] * other.m[3] + m[13] * other.m[7] + m[14] * other.m[11] + m[15] * other.m[15]) % k;
            memcpy(m, buffer, sizeof(unsigned int) * n * n);
            return;
        }

        for (unsigned int r = 0; r < n; ++r) {
            for (unsigned int c = 0; c < n; ++c) {
                unsigned int sum = 0;
                for (unsigned int i = 0; i < n; ++i) {
                    sum += m[r * n + i] * other.m[i * n + c];
                }
                buffer[r * n + c] = sum % k;
            }
        }
        memcpy(m, buffer, sizeof(unsigned int) * n * n);
    }

    void power(mpz_class e, Matrix &result) {
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
        for (unsigned int r = 0; r < n; ++r) {
            for (unsigned int c = 0; c < n; ++c) {
                printf("%8d", at(r, c));
                if (c == n - 1) {
                    printf("\n");
                } else {
                    printf("\t");
                }
            }
        }
    }

 private:
    const unsigned int n, k;
    unsigned int m[20 * 20];
    static unsigned int buffer[20 * 20];
};

unsigned int Matrix::buffer[20 * 20];

}

#endif
