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
#include "Mutex.h"

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
    class MatrixIterator {
     public:
        MatrixIterator()
            :done(false) {
        }

        virtual void next() = 0;

        const Matrix &matrix() {
            return m;
        }

        bool is_done() const {
            return done;
        }

     protected:
        Matrix<N, K> m;
        bool done;
    };

    class GeneralIterator : public MatrixIterator {
     public:
        GeneralIterator() : MatrixIterator() {
            MatrixIterator::m.zero();
        }

        virtual void next() {
            if (!MatrixIterator::m.next()) {
                MatrixIterator::done = true;
            }
        }
    };

    class RandomIterator : public MatrixIterator {
     public:
        RandomIterator(unsigned int size) : MatrixIterator(),
            size(size) {
        }

        virtual void next() {
            if (MatrixIterator::done) {
                return;
            }

            MatrixIterator::m.randomize();
            --size;
            if (size == 0) {
                MatrixIterator::done = true;
            }
        }

     protected:
        unsigned int size;
    };


    class NonsingularIterator : public MatrixIterator {
        typedef vector<unsigned int> Row;

     public:
        NonsingularIterator(bool non_permutating) : MatrixIterator(),
            non_permutating(non_permutating) {
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
            MatrixIterator::done = false;
        }

        virtual void next() {
            if (!MatrixIterator::done) {
                next_matrix(N - 1);
            }
        }

     protected:
        void next_matrix(unsigned int row_number) {
            if (MatrixIterator::done) {
                return;
            }
            if (row_number > 0 && state[row_number - 1] == -1) {
                next_matrix(row_number - 1);
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
                memcpy(&MatrixIterator::m.m[row_number * N], row.data(), sizeof(unsigned int) * N);
                if (row_number < N - 1) {
                    rule_out_rows(row_number);
                }
                return;
            }

            if (row_number == 0) {
                MatrixIterator::done = true;
            } else {
                state[row_number] = -1;
                next_matrix(row_number - 1);
                next_matrix(row_number);
            }
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

        vector<unsigned int> coefficients;
        vector<unsigned int> base_rows;
        vector<Row> all_rows;
        vector<bool> acceptable_rows;
        vector<vector<unsigned int>> ruled_out_rows;
        vector<unsigned int> base_ruled_out_rows;
        int state[N];
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

    inline void mult(const Matrix<N, K> &other, unsigned int **external_buffer) {
        unsigned int *buffer_used = buffer;
        if (external_buffer) {
            buffer_used = *external_buffer;
        }
        for (unsigned int r = 0; r < N; ++r) {
            for (unsigned int c = 0; c < N; ++c) {
                unsigned int sum = 0;
                for (unsigned int i = 0; i < N; ++i) {
                    sum += at(r, i) * other.at(i, c);
                }
                buffer_used[r * N + c] = sum % K;
            }
        }
        if (external_buffer) {
            unsigned int *tmp = *external_buffer;
            *external_buffer = m;
            m = tmp;
        } else {
            unsigned int *tmp = buffer;
            buffer = m;
            m = tmp;
        }
#ifdef DEBUG
        ++num_multiplications;
#endif
    }

    inline void square(unsigned int **external_buffer) {
        unsigned int *buffer_used = buffer;
        if (external_buffer) {
            buffer_used = *external_buffer;
        }
        for (unsigned int r = 0; r < N; ++r) {
            for (unsigned int c = 0; c < N; ++c) {
                unsigned int sum = 0;
                for (unsigned int i = 0; i < N; ++i) {
                    sum += at(r, i) * at(i, c);
                }
                buffer_used[r * N + c] = sum % K;
            }
        }
        if (external_buffer) {
            unsigned int *tmp = *external_buffer;
            *external_buffer = m;
            m = tmp;
        } else {
            unsigned int *tmp = buffer;
            buffer = m;
            m = tmp;
        }
#ifdef DEBUG
        ++num_multiplications;
#endif
    }

    void power(const mpz_class &e, Matrix<N, K> &result, unsigned int **external_buffer) const {
        if (e == 0) {
            result.make_identity();
            return;
        } else if (e == 1) {
            result = *this;
            return;
        }

        power(e / 2, result, external_buffer);
        result.square(external_buffer);
        if (e % 2 == 1) {
            result.mult(*this, external_buffer);
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

    Cycle get_cycle(const mpz_class &size, set<mpz_class> &seen, unsigned int **external_buffer) const {
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
            tmp.mult(*this, external_buffer);
        }
        cycle.push_back(tmp.get_id());

        return cycle;
    }

    static unsigned int search_for_prime_exponent(const Matrix<N, K> &base, const Factor &factor,
        unsigned int **external_buffer) {
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
            base.power(exp, tmp, external_buffer);
            if (tmp.is_identity()) {
                return result;
            }
        }
        return result;
    }

    mpz_class get_minimal_exponent(const mpz_class &phi, const Factors &factors, Factors *left_factors,
        unsigned int **external_buffer) {
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
            power(exponent, tmp, external_buffer);
            unsigned int exp = search_for_prime_exponent(tmp, factors[i], external_buffer);
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
        bool check_all = false;
        if (phi < 50000000) {
            check_all = true;
        }
        unsigned int z = 0;
        unsigned int s = 0;
        unsigned int c = 0;
        mpz_class max_c = 2000000, ctmp;
        if (check_all) {
            max_c = phi;
        }
        unsigned int last_timestamp = time(NULL);

        bool goon = true;
        Mutex lock_iterator;
        Mutex lock_result;
        mpz_class result = 1;
        bool result_changed = false;

        #pragma omp parallel firstprivate(factors)
        {
            MatrixIterator *it;
            if (check_all) {
                //it = new GeneralIterator();
                it = new NonsingularIterator(true);
            } else {
                it = new RandomIterator(1000000000);
            }
            mpz_class local_result = 1;
            Matrix<N, K> tmp, tmp2;
            Matrix<N, K> m;
            Matrix<N, K> base;
            unsigned int *buffer = new unsigned int [N * N];
            while (goon) {
                //lock_iterator.lock();
                goon = !it->is_done();
                //lock_iterator.unlock();
                if (!goon) {
                    break;
                }
                //lock_iterator.lock();
                it->next();
                m = it->matrix();
                //lock_iterator.unlock();

                if (c > max_c) {
                    break;
                }

                #pragma omp master
                {
                    if (time(NULL) >= last_timestamp + 5) {
                        ctmp = 10000 * mpz_class(c) / max_c;
                        printf("%.f%% c:%d, z:%d, s:%d\n", (float) ctmp.get_ui() / 100, c, z, s);
                        last_timestamp = time(NULL);
                    }

                    lock_result.lock();
                    if (result_changed) {
                        printf("%s\n", result.get_str().c_str());
                        result_changed = false;
                    }
                    lock_result.unlock();
                }

                mpz_class e = 1;
                m.power(phi, tmp, &buffer);
                if (tmp.is_zero()) {
                    #pragma omp atomic
                    ++z;
                    continue;
                } else if (!tmp.is_identity()) {
                    #pragma omp atomic
                    ++s;
                    // singular matrix, ignore
                    continue;
                }

                m.power(local_result, tmp, &buffer);
                if (tmp.is_identity()) {
                    #pragma omp atomic
                    ++c;
                    continue;
                }

                /*for (unsigned int i = 0; i < factors.size(); ++i) {
                    printf("%s^%d, ", factors[i].prime.get_str().c_str(), factors[i].exponent);
                }
                printf("\n");*/
                mpz_class rest;
                rest = phi / local_result;
                e = tmp.get_minimal_exponent(rest, factors, &factors, &buffer);
                /*for (unsigned int i = 0; i < factors.size(); ++i) {
                    printf("%s^%d, ", factors[i].prime.get_str().c_str(), factors[i].exponent);
                }
                printf("\n");*/
                #pragma omp atomic
                ++c;

                if (e > 1) {
                    lock_result.lock();
                    local_result *= e;
                    mpz_class old_result = result;
                    mpz_lcm(result.get_mpz_t(), result.get_mpz_t(), local_result.get_mpz_t());
                    if (result != old_result) {
                        result_changed = true;
                    }
                    lock_result.unlock();
                }
                /*unsigned int prime_factors = 0;
                for (unsigned int i = 0; i < factors.size(); ++i) {
                    prime_factors += factors[i].exponent;
                }
                printf("rest: %d\n", prime_factors);
                */
            }

            delete[] buffer;
        }

        unsigned int confidence;
        if (check_all) {
            confidence = 0;
        } else {
            confidence = c;
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

            m.power(phi, tmp, nullptr);
            if (tmp.is_zero()) {
                ++z;
                continue;
            } else if (!tmp.is_identity()) {
                ++s;
                // singular matrix, ignore
                continue;
            }

            ++c;
            mpz_class e = m.get_minimal_exponent(phi, factors, nullptr, nullptr);
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
                cycles.push_back(m.get_cycle(e, seen, nullptr));
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
