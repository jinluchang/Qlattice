#pragma once

inline int omp_get_max_threads()
{
  return 1;
}

inline int omp_get_num_threads()
{
  return 1;
}

inline int omp_get_thread_num()
{
  return 0;
}
