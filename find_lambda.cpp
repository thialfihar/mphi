#include <cstdio>
#include <cstring>

#include <algorithm>
#include <vector>
#include <gmpxx.h>

#include "Matrix.h"

void test_stuff() {
    Matrix<3> m1(7);
    Matrix<3> m2(7);
    Matrix<3> result(7);

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
            mpz_class lambda;
            if (n == 2) {
                lambda = Matrix<2>::find_probable_lambda(k);
            } else if (n == 3) {
                lambda = Matrix<3>::find_probable_lambda(k);
            } else if (n == 4) {
                lambda = Matrix<4>::find_probable_lambda(k);
            } else if (n == 5) {
                lambda = Matrix<5>::find_probable_lambda(k);
            } else if (n == 6) {
                lambda = Matrix<6>::find_probable_lambda(k);
            } else if (n == 7) {
                lambda = Matrix<7>::find_probable_lambda(k);
            } else if (n == 8) {
                lambda = Matrix<8>::find_probable_lambda(k);
            }
            printf("lambda(n=%d, k=%d) = %s\n", n, k, lambda.get_str().c_str());
            fflush(stdout);
        }
    }
    return 0;
}
