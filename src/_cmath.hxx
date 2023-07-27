#pragma once
#include <type_traits>
#include <cmath>
#include <random>

using std::is_floating_point;
using std::uniform_int_distribution;
using std::ceil;
using std::sqrt;




#pragma region METHODS
#pragma region CEIL DIV
/**
 * Calculate ceil(x/y).
 * @param x the dividend
 * @param y the divisor
 * @returns ceil(x/y)
 */
template <class T>
inline T ceilDiv(T x, T y) {
  // For kernel launch calculation.
  if (is_floating_point<T>::value) return ceil(x/y);
  return (x + y-1) / y;
}
#pragma endregion




#pragma region SGN
/**
 * Get the sign of a value.
 * @param x the value
 * @returns -1 if x<0, 0 if x==0, 1 if x>0
 */
template <typename T>
inline int sgn(T x) {
  return (T() < x) - (x < T());
}
// - https://stackoverflow.com/a/4609795/1413259
#pragma endregion




#pragma region POW2
/**
 * Check if a value is a power of 2.
 * @param x the value
 * @returns true if x is a power of 2
 */
template <class T>
constexpr bool isPow2(T x) noexcept {
  return !(x & (x-1));
}


/**
 * Get the previous power of 2 of a value.
 * @param x the value
 * @returns previous power of 2 of x
 */
template <class T>
constexpr T prevPow2(T x) noexcept {
  return 1 << T(log2(x));
}


/**
 * Get the next power of 2 of a value.
 * @param x the value
 * @returns next power of 2 of x
 */
template <class T>
constexpr T nextPow2(T x) noexcept {
  return 1 << T(ceil(log2(x)));
}
#pragma endregion




#pragma region PRIME
/**
 * Examine if a value is prime.
 * @param x the value
 * @returns true if x is prime
 */
template <class T>
inline bool isPrime(T x) {
  // 1. 2, 3 are prime
  if (x<=3) return x>1;
  // 2. Multiples of 2, 3 not prime
  if (x % 2==0 || x % 3==0) return false;
  // 3. Factor of 6k-1 or 6k+1 => not prime
  for (T i=6, I=T(sqrt(x))+1; i<=I; i+=6)
    if (x % (i-1)==0 || x % (i+1)==0) return false;
  return true;
}


/**
 * Obtain the next prime of a value.
 * @param x the value
 * @returns next prime of x
 */
template <class T>
inline T nextPrime(T x) {
  while (true)
    if (isPrime(++x)) return x;
}


/**
 * Obtain a random prime number in a range.
 * @param begin beginning of the range
 * @param end end of the range
 * @param rnd random number generator
 */
template <class T, class R>
inline T randomPrime(T begin, T end, R& rnd) {
  uniform_int_distribution<T> dis(begin, end);
  for (int i=128; i>0; --i) {
    T a = dis(rnd);
    if (isPrime(a)) return a;
  }
  return end-1;
}
#pragma endregion
#pragma endregion
