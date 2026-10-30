#pragma once

#ifndef _OPENMP
inline int omp_get_num_threads() { return 1; }
inline void omp_set_num_threads(int) {}
inline void omp_set_max_active_levels(int) {}
inline void omp_set_nested(int) {}
#endif
