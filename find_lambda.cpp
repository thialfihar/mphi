#include <cstdio>
#include <cstring>

#include "Matrix.h"

#define N 9
#define K 3

int main(int argc, char *argv[]) {
    srand(time(nullptr));
    //test_stuff(); return 0;

    mpz_class lambda = Matrix<N, K>::find_probable_lambda();
    printf("lambda(n=%d, k=%d) = %s\n", N, K, lambda.get_str().c_str());
    return 0;
}
