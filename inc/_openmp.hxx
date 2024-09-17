#pragma once
#include <omp.h>




#pragma region BELONGS
/**
 * Check if work belongs to current thread.
 * @param key work key
 * @param thread current thread
 * @param THREADS available threads
 * @returns true if work belongs to current thread
 */
template <class K>
inline bool belongsOmp(K key, int thread, int THREADS) {
  const K CHUNK_SIZE = 1024;
  K chunk = key / CHUNK_SIZE;
  return chunk % THREADS == K(thread);
}


/**
 * Check if work belongs to current thread.
 * @param key work key
 * @returns true if work belongs to current thread
 */
template <class K>
inline bool belongsOmp(K key) {
  int thread  = omp_get_thread_num();
  int THREADS = omp_get_num_threads();
  return belongsOmp(key, thread, THREADS);
}
#pragma endregion
