#ifndef __MUTEX_H__
#define __MUTEX_H__

#ifdef _OPENMP

# include <omp.h>
class Mutex {
 public:
    Mutex() {
        omp_init_lock(&omp_lock);
        locked = false;
    }

    Mutex(const Mutex&) : Mutex() {
    }

    ~Mutex() {
        unlock();
        omp_destroy_lock(&omp_lock);
    }

    void lock() {
        omp_set_lock(&omp_lock);
        locked = true;
    }

    void unlock() {
        if (!locked) {
            return;
        }
        locked = false;
        omp_unset_lock(&omp_lock);
    }

 protected:
    omp_lock_t omp_lock;
    bool locked;
};

#else

class Mutex {
 public:
    void lock() {}
    void unlock() {}
};

#endif


#endif //  __MUTEX_H__
