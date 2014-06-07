#ifndef __MATRIX_H__
#define __MATRIX_H__
#include <assert.h>
#include <cstdio>
#include <cstring>
#include <ctime>

#include <algorithm>
#include <iterator>
#include <map>
#include <set>
#include <vector>
#include <gmpxx.h>

#include "functions.h"

//#define DEBUG

using std::set;
using std::map;

typedef struct {
    mpz_class lambda;
    unsigned int confidence;
    mpz_class phi;
} LambdaResult;

typedef std::map<mpz_class, mpz_class> ExponentMap;
typedef vector<mpz_class> Cycle;

typedef struct {
    ExponentMap map;
    ExponentMap permutation_map;
    unsigned int confidence;
    mpz_class phi;
    vector<Cycle> cycles;
} ExponentMapResult;

bool cmp_factors(const Factor &a, const Factor &b) {
    mpz_class t1, t2;
    mpz_pow_ui(t1.get_mpz_t(), a.prime.get_mpz_t(), a.exponent);
    mpz_pow_ui(t2.get_mpz_t(), b.prime.get_mpz_t(), b.exponent);
    return t1 > t2;
}

template <unsigned int N, unsigned int K>
class Matrix {
 public:
    class GeneralIterator : public std::iterator<std::input_iterator_tag, Matrix<N, K>> {
     public:
        GeneralIterator() {
            m.zero();
            done = false;
        }

        GeneralIterator(const GeneralIterator &o)
            :m(o.m), done(o.done) {
        }

        GeneralIterator &operator++() {
            if (!m.next()) {
                done = true;
            }
            return *this;
        }

        bool operator==(const GeneralIterator &o) {
            return m == o.m;
        }

        bool operator!=(const GeneralIterator &o) {
            return !(*this == o);
        }

        const Matrix &operator*() {
            return m;
        }

        const Matrix *operator&() {
            return &m;
        }

        bool is_done() const {
            return done;
        }

     private:
        Matrix<N, K> m;
        bool done;
    };

    class NonsingularIterator : public std::iterator<std::input_iterator_tag, Matrix<N, K>> {
        typedef vector<unsigned int> Row;

     public:
        NonsingularIterator() {
            Row row(N, 0);
            ruled_out_rows = vector<vector<unsigned int>>(N);
            unsigned int max_rows = 1;
            for (unsigned int i = 0; i < N; ++i) {
                max_rows *= K;
                // reserve more and more for each level
                ruled_out_rows[i].reserve(max_rows);
            }

            for (unsigned int i = 1; i < K; ++i) {
                if (gcd(i, K) == 1) {
                    coefficients.push_back(i);
                }
            }

            all_rows.push_back(row);
            while (true) {
                int i = N - 1;
                bool found_one = false;
                while (i >= 0) {
                    row[i] = (row[i] + 1) % K;
                    if (row[i] != 0) {
                        found_one = true;
                        break;
                    }
                    --i;
                }
                if (!found_one) {
                    break;
                }
                all_rows.push_back(row);
            }

            acceptable_rows = vector<bool>(all_rows.size(), false);
            for (unsigned int i = 0; i < all_rows.size(); ++i) {
                unsigned int g = K;
                for (unsigned int j = 0; j < row.size(); ++j) {
                    g = gcd(g, all_rows[i][j]);
                    if (g == 1) {
                        break;
                    }
                }
                if (g == 1) {
                    base_rows.push_back(i);
                    acceptable_rows[i] = true;
                } else {
                    base_ruled_out_rows.push_back(i);
                }
            }

            printf("rows: %d, usable: %d\n", (unsigned int) all_rows.size(),
                (unsigned int) base_rows.size());
            for (unsigned int i = 0; i < N; ++i) {
                state[i] = -1;
            }
            done = false;
        }

        NonsingularIterator(const NonsingularIterator &o)
            :m(o.m), done(o.done) {
        }

        void clear_ruled_out_rows(unsigned int row_number) {
            vector<unsigned int> &relevant = ruled_out_rows[row_number];
            for (unsigned int i = 0; i < relevant.size(); ++i) {
                acceptable_rows[relevant[i]] = true;
            }
            relevant.clear();
        }

        void rule_out_rows(unsigned int row_number) {
            // this contains the 0 row and for composite K all non relative prime multiples of
            // rows
            for (unsigned int j = 0; j < base_ruled_out_rows.size(); ++j) {
                const Row &row = all_rows[base_ruled_out_rows[j]];
                for (unsigned int c = 0; c < coefficients.size(); ++c) {
                    unsigned int tmp = 0;
                    for (unsigned int ii = 0; ii < N; ++ii) {
                        tmp *= K;
                        tmp += (row[ii] + all_rows[base_rows[state[row_number]]][ii] * coefficients[c]) % K;

                    }
                    if (acceptable_rows[tmp]) {
                        acceptable_rows[tmp] = false;
                        ruled_out_rows[row_number].push_back(tmp);
                    }
                }
            }

            for (unsigned int i = 0; i < row_number; ++i) {
                for (unsigned int j = 0; j < ruled_out_rows[i].size(); ++j) {
                    const Row &row = all_rows[ruled_out_rows[i][j]];
                    for (unsigned int c = 0; c < coefficients.size(); ++c) {
                        unsigned int tmp = 0;
                        for (unsigned int ii = 0; ii < N; ++ii) {
                            tmp *= K;
                            tmp += (row[ii] + all_rows[base_rows[state[row_number]]][ii] * coefficients[c]) % K;

                        }
                        if (acceptable_rows[tmp]) {
                            acceptable_rows[tmp] = false;
                            ruled_out_rows[row_number].push_back(tmp);
                        }
                    }
                }
            }
        }

        void next(unsigned int row_number) {
            if (done) {
                return;
            }
            if (row_number > 0 && state[row_number - 1] == -1) {
                next(row_number - 1);
            }

            clear_ruled_out_rows(row_number);

            if (non_permutating && state[row_number] == -1 && row_number > 0) {
                state[row_number] = state[row_number - 1];
                //printf("new: %d:%d %d:%d\n", row_number, state[row_number], row_number - 1, state[row_number - 1]);
            }
            for (int i = state[row_number] + 1; i < (int) base_rows.size(); ++i) {
                //printf("hmm %d\n", i);
                const unsigned int id = base_rows[i];
                if (!acceptable_rows[id]) {
                    continue;
                }
                const Row &row = all_rows[id];

                state[row_number] = i;
                memcpy(&m.m[row_number * N], row.data(), sizeof(unsigned int) * N);
                if (row_number < N - 1) {
                    rule_out_rows(row_number);
                }
                return;
            }

            if (row_number == 0) {
                done = true;
            } else {
                state[row_number] = -1;
                next(row_number - 1);
                next(row_number);
            }
        }

        NonsingularIterator &operator++() {
            if (!done) {
                next(N - 1);
            }
            return *this;
        }

        bool operator==(const NonsingularIterator &o) {
            return m == o.m;
        }

        bool operator!=(const NonsingularIterator &o) {
            return !(*this == o);
        }

        const Matrix &operator*() {
            return m;
        }

        const Matrix *operator&() {
            return &m;
        }

        bool is_done() const {
            return done;
        }

        void set_non_permutating() {
            non_permutating = true;
        }

     private:
        Matrix<N, K> m;
        vector<unsigned int> coefficients;
        vector<unsigned int> base_rows;
        vector<Row> all_rows;
        vector<bool> acceptable_rows;
        vector<vector<unsigned int>> ruled_out_rows;
        vector<unsigned int> base_ruled_out_rows;
        int state[N];
        bool done;
        bool non_permutating = false;
    };

    inline Matrix() {
        m = new unsigned int[N * N];
    }

    inline Matrix(const Matrix<N, K> &other)
        :Matrix() {
        *this = other;
    }

    inline ~Matrix() {
        if (m) {
            delete[] m;
        }
    }

    inline void operator==(const Matrix<N, K> &other) {
        for (unsigned int i = 0; i < N * N; ++i) {
            if (m[i] != other.m) {
                return false;
            }
        }
        return true;
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

    inline bool is_permutation_matrix() const {
        for (unsigned int r = 0; r < N; ++r) {
            unsigned int ones = 0;
            unsigned int zeros = 0;
            for (unsigned int c = 0; c < N; ++c) {
                unsigned int value = at(r, c);
                if (value == 1) {
                    ++ones;
                } else if (value == 0) {
                    ++zeros;
                }
            }

            if (ones != 1 || zeros != N - 1) {
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

    mpz_class get_id() const {
        mpz_class id = 0;
        for (unsigned int i = 0; i < N * N; ++i) {
            unsigned int j = N * N - 1 - i;
            id *= K;
            id += m[j];
        }
        return id;
    }

    Cycle get_cycle(const mpz_class &size, set<mpz_class> &seen) const {
        Cycle cycle;
        assert(size.fits_uint_p());
        if (size != 0) {
            cycle.reserve(size.get_ui());
        }
        Matrix<N, K> tmp = *this;
        while (!tmp.is_identity()) {
            mpz_class id = tmp.get_id();
            seen.insert(id);
            cycle.push_back(id);
            tmp.mult(*this);
        }
        cycle.push_back(tmp.get_id());

        return cycle;
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
        Factors new_factors;
        if (left_factors) {
            new_factors = factors;
        }

        for (unsigned int i = 0; i < factors.size(); ++i) {
            for (unsigned int j = 0; j < factors[i].exponent; ++j) {
                exponent /= factors[i].prime;
            }
            power(exponent, tmp);
            unsigned int exp = search_for_prime_exponent(tmp, factors[i]);
            //printf("%s : %d : %d\n", factors[i].prime.get_str().c_str(), factors[i].exponent, exp);
            if (left_factors) {
                new_factors[i].exponent -= exp;
            }
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
        unsigned int z = 0;
        unsigned int s = 0;
        mpz_class c, max_c = 2000000, ctmp;
        if (check_all) {
            max_c = phi;
        }
        Matrix<N, K> tmp, tmp2;
        unsigned int last_timestamp = time(NULL);
        const Matrix<N, K> *m;
        Matrix<N, K> base;
        NonsingularIterator it;
        it.set_non_permutating();
        while (true) {
            if (check_all) {
                ++it;
                if (it.is_done()) {
                    break;
                }
                m = &it;
            } else {
                if (c >= max_c) {
                    break;
                }
                base.randomize();
                m = &base;
            }
            if (time(NULL) >= last_timestamp + 5) {
                ctmp = 10000 * c / max_c;
                printf("%.f%% c:%s, z:%d, s:%d\n", (float) ctmp.get_ui() / 100, c.get_str().c_str(), z, s);
                last_timestamp = time(NULL);
            }
            /*++c;
            if (c >= 20000000) {
                m->print();
            m->power(phi, tmp);
            tmp.print();
                break;
            }
            continue;*/

            //printf("%d %d\n", c, z);
            //m.print();
            mpz_class e = 1;
            m->power(phi, tmp);
            if (tmp.is_zero()) {
                ++z;
                continue;
            } else if (!tmp.is_identity()) {
                ++s;
                // singular matrix, ignore
                continue;
            }

            m->power(result, tmp);
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
        return {result, confidence, phi};
    }

    static ExponentMapResult get_exponent_map() {
        mpz_class phi = phi_n(N, K);
        Factors factors = factorize(phi);
        // go backwards, so large factors are removed early on
        sort(factors.begin(), factors.end(), cmp_factors);
        ExponentMap result;
        ExponentMap permutation_map;
        vector<Cycle> cycles;
        set<mpz_class> seen;

        mpz_class num_matrices = 1;
        for (unsigned int i = 0; i < N * N; ++i) {
            num_matrices *= K;
        }
        bool check_all = false;
        if (phi < 50000000) {
            check_all = true;
        }
        Matrix<N, K> m;
        m.zero();
        unsigned int z = 0;
        unsigned int s = 0;
        mpz_class c = 0, max_c = 200000, ctmp;
        if (check_all) {
            max_c = phi;
        }
        cycles.reserve(max_c.get_ui());
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

            m.power(phi, tmp);
            if (tmp.is_zero()) {
                ++z;
                continue;
            } else if (!tmp.is_identity()) {
                ++s;
                // singular matrix, ignore
                continue;
            }

            ++c;
            mpz_class e = m.get_minimal_exponent(phi, factors, nullptr);
            auto it = result.find(e);
            if (it == result.end()) {
                result[e] = 1;
                if (m.is_permutation_matrix()) {
                    permutation_map[e] = 1;
                }
            } else {
                ++it->second;
                if (m.is_permutation_matrix()) {
                    ++permutation_map[e];
                }
            }

            mpz_class id = m.get_id();
            if (seen.find(id) == seen.end()) {
                cycles.push_back(m.get_cycle(e, seen));
            }
        }

        unsigned int confidence;
        if (check_all) {
            confidence = 0;
        } else {
            confidence = c.get_ui();
        }
        return {result, permutation_map, confidence, phi, cycles};
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
