#include <cstdio>
#include <cstring>
#include <ctime>

#include "Matrix.h"
#include "GeneralMatrix.h"

#define MAX_MULTIPLICATIONS 1000000
#define MAX_POWERS 100000
#define N 10
#define K 10

float timediff(clock_t t1, clock_t t2) {
    return (float) (t2 - t1) / CLOCKS_PER_SEC;
}

int main(int argc, char *argv[]) {
    srand(2107);

    Matrix<N, K> m1, m2, result;
    m1.randomize();
    m2.randomize();
    clock_t t1 = clock();
    clock_t noise_t = 0;
    for (unsigned int i = 0; i < MAX_MULTIPLICATIONS; ++i) {
        m1.next();
        m1.mult(m2, result);
    }
    clock_t t2 = clock();
    t1 += noise_t;
    float dt = timediff(t1, t2);
    printf("Matrix<%d>(%d) mult x %d: %.3f (%.3f m/s)\n", N, K, MAX_MULTIPLICATIONS, dt, MAX_MULTIPLICATIONS / dt);

    gm::Matrix gm1(N, K), gm2(N, K), gresult(N, K);
    gm1.randomize();
    gm2.randomize();
    t1 = clock();
    noise_t = 0;
    for (unsigned int i = 0; i < MAX_MULTIPLICATIONS; ++i) {
        gm1.next();
        gm1.mult(gm2, gresult);
    }
    t2 = clock();
    t1 += noise_t;
    dt = timediff(t1, t2);
    printf("Matrix(%d, %d) mult x %d: %.3f (%.3f m/s)\n", N, K, MAX_MULTIPLICATIONS, dt, MAX_MULTIPLICATIONS / dt);

    printf("\n");

    t1 = clock();
    noise_t = 0;
    for (unsigned int i = 0; i < MAX_MULTIPLICATIONS; ++i) {
        m1.next();
        m1.mult(m2);
    }
    t2 = clock();
    t1 += noise_t;
    dt = timediff(t1, t2);
    printf("Matrix<%d>(%d) mult in place x %d: %.3f (%.3f m/s)\n", N, K, MAX_MULTIPLICATIONS, dt, MAX_MULTIPLICATIONS / dt);

    t1 = clock();
    noise_t = 0;
    for (unsigned int i = 0; i < MAX_MULTIPLICATIONS; ++i) {
        gm1.next();
        gm1.mult(gm2);
    }
    t2 = clock();
    t1 += noise_t;
    dt = timediff(t1, t2);
    printf("Matrix(%d, %d) mult in place x %d: %.3f (%.3f m/s)\n", N, K, MAX_MULTIPLICATIONS, dt, MAX_MULTIPLICATIONS / dt);

    printf("\n");

    m1.randomize();
    t1 = clock();
    noise_t = 0;
    unsigned int e = 174762; // binary: 101010101010101010
    for (unsigned int i = 0; i < MAX_POWERS; ++i) {
        m1.next();
        m1.power(e, result);
    }
    t2 = clock();
    t1 += noise_t;
    dt = timediff(t1, t2);
    printf("Matrix<%d>(%d) power x %d: %.3f (%.3f p/s)\n", N, K, MAX_POWERS, dt, MAX_POWERS / dt);

    gm1.randomize();
    t1 = clock();
    noise_t = 0;
    for (unsigned int i = 0; i < MAX_POWERS; ++i) {
        gm1.next();
        gm1.power(e, gresult);
    }
    t2 = clock();
    t1 += noise_t;
    dt = timediff(t1, t2);
    printf("Matrix(%d, %d) power x %d: %.3f (%.3f p/s)\n", N, K, MAX_POWERS, dt, MAX_POWERS / dt);


    return 0;
}
