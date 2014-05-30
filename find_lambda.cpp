#include <cstdio>
#include <cstring>
#include <boost/preprocessor/repetition/repeat.hpp>

#include "Matrix.h"

#define MIN_N 2
#define NUM_N 7
#define MIN_K 2
#define NUM_K 7

#define NEXT_N(z, N, K) else if (n == MIN_N + N && k == MIN_K + K) { result = Matrix<MIN_N + N, MIN_K + K>::find_probable_lambda(); }
#define NEXT_K(z, K, data) BOOST_PP_REPEAT_2(NUM_N, NEXT_N, K)

int main(int argc, char *argv[]) {
    srand(time(nullptr));
    //test_stuff(); return 0;

    for (unsigned int k = MIN_K; k < MIN_K + NUM_K; ++k) {
        for (unsigned int n = MIN_N; n < MIN_N + NUM_N; ++n) {
            LambdaResult result;
            if (false) {
                // ...
            }
            BOOST_PP_REPEAT_1(NUM_K, NEXT_K, 0);
            printf("lambda(n=%d, k=%d) = %s [%d:%s]\n",
                n, k, result.lambda.get_str().c_str(), result.confidence, result.phi.get_str().c_str());
        }
    }
    return 0;
}
