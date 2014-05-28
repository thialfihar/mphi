#include <cstdio>
#include <cstring>

#include <algorithm>
#include <vector>
#include <gmpxx.h>

#include "Matrix.h"

int main(int argc, char *argv[]) {
    srand(time(nullptr));
    //test_stuff(); return 0;

    mpz_class lambda = Matrix<6, 6>::find_probable_lambda();
    printf("lambda(n=%d, k=%d) = %s\n", 6, 6, lambda.get_str().c_str());
    return 0;
}
