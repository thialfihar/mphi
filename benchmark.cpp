#include <cstdio>
#include <cstring>
#include <ctime>

#include "Matrix.h"
#include "GeneralMatrix.h"

#define MAX_MULTIPLICATIONS 1000000
#define N 6
#define K 65

float timediff(clock_t t1, clock_t t2) {
    return (float) (t2 - t1) / CLOCKS_PER_SEC;
}

int main(int argc, char *argv[]) {
    srand(2107);

    clock_t t1 = clock();
    Matrix<N> m1(K), m2(K), result(K);
    m1.randomize();
    m2.randomize();
    clock_t noise_t = 0;
    for (unsigned int i = 0; i < MAX_MULTIPLICATIONS / 10; ++i) {
        /*clock_t tmp = clock();
        m1.randomize();
        m2.randomize();
        noise_t += clock() - tmp;*/
        m1.mult(m2, result);
        m1.mult(m2, result);
        m1.mult(m2, result);
        m1.mult(m2, result);
        m1.mult(m2, result);
        m1.mult(m2, result);
        m1.mult(m2, result);
        m1.mult(m2, result);
        m1.mult(m2, result);
        m1.mult(m2, result);
    }
    clock_t t2 = clock();
    t1 += noise_t;
    float dt = timediff(t1, t2);
    printf("Matrix<%d>(%d) mult x %d: %.3f (%.3f m/s)\n", N, K, MAX_MULTIPLICATIONS, dt, MAX_MULTIPLICATIONS / dt);

    t1 = clock();
    gm::Matrix gm1(N, K), gm2(N, K), gresult(N, K);
    gm1.randomize();
    gm2.randomize();
    noise_t = 0;
    for (unsigned int i = 0; i < MAX_MULTIPLICATIONS / 10; ++i) {
        /*clock_t tmp = clock();
        gm1.randomize();
        gm2.randomize();
        noise_t += clock() - tmp;*/
        gm1.mult(gm2, gresult);
        gm1.mult(gm2, gresult);
        gm1.mult(gm2, gresult);
        gm1.mult(gm2, gresult);
        gm1.mult(gm2, gresult);
        gm1.mult(gm2, gresult);
        gm1.mult(gm2, gresult);
        gm1.mult(gm2, gresult);
        gm1.mult(gm2, gresult);
        gm1.mult(gm2, gresult);
    }
    t2 = clock();
    t1 += noise_t;
    dt = timediff(t1, t2);
    printf("Matrix(%d, %d) mult x %d: %.3f (%.3f m/s)\n", N, K, MAX_MULTIPLICATIONS, dt, MAX_MULTIPLICATIONS / dt);


    return 0;
}
