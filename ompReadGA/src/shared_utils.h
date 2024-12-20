#ifndef SHARED_UTILS_H
#define SHARED_UTILS_H

#ifdef _OPENMP
#include <omp.h>
#endif

inline unsigned int use_threads(int n_threads = 1) {
#ifdef _OPENMP
  if (n_threads <= 1) return 1;
  if (n_threads > omp_get_max_threads()) {
    return static_cast<unsigned int>(omp_get_max_threads());
  }
  return static_cast<unsigned int>(n_threads);
#else
  return 1;
#endif
}

#endif // SHARED_UTILS_H
