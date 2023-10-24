#pragma once
#include <utility>
#include <type_traits>
#include <iterator>
#include <array>
#include <vector>
#include <ostream>
#include <iostream>
#include <chrono>
#include <ctime>

using std::pair;
using std::array;
using std::vector;
using std::ostream;
using std::is_fundamental;
using std::iterator_traits;
using std::time_t;
using std::tm;
using std::localtime;
using std::chrono::time_point;
using std::chrono::system_clock;
using std::cout;




#pragma region WRITE
/**
 * Write values to a stream.
 * @param a the stream
 * @param ib begin iterator of values
 * @param ie end iterator of values
 */
template <class I>
inline void write_values(ostream& a, I ib, I ie) {
  using T = typename iterator_traits<I>::value_type;
  if (is_fundamental<T>::value) {
    a << "{";
    for (; ib < ie; ++ib)
      a << " " << *ib;
    a << " }";
  }
  else {
    a << "{\n";
    for (; ib < ie; ++ib)
      a << "  " << *ib << "\n";
    a << "}";
  }
}


/**
 * Write values to a stream.
 * @param a the stream
 * @param x the container of values
 */
template <class J>
inline void writeValues(ostream& a, const J& x) {
  write_values(a, x.begin(), x.end());
}


/**
 * Write a pair to a stream.
 * @param a the stream
 * @param x the pair
 */
template <class K, class V>
inline void write(ostream& a, const pair<K, V>& x) {
  a << x.first << ": " << x.second;
}


/**
 * Write an array to a stream.
 * @param a the stream
 * @param x the array
 */
template <class T, size_t N>
inline void write(ostream& a, const array<T, N>& x) {
  writeValues(a, x);
}


/**
 * Write a vector to a stream.
 * @param a the stream
 * @param x the vector
 */
template <class T>
inline void write(ostream& a, const vector<T>& x) {
  writeValues(a, x);
}


/**
 * Write a pair to a stream.
 * @param a the stream
 * @param x the pair
 */
template <class K, class V>
inline ostream& operator<<(ostream& a, const pair<K, V>& x) {
  write(a, x); return a;
}


/**
 * Write an array to a stream.
 * @param a the stream
 * @param x the array
 */
template <class T, size_t N>
inline ostream& operator<<(ostream& a, const array<T, N>& x) {
  write(a, x); return a;
}


/**
 * Write a vector to a stream.
 * @param a the stream
 * @param x the vector
 * @returns the stream
 */
template <class T>
inline ostream& operator<<(ostream& a, const vector<T>& x) {
  write(a, x); return a;
}
#pragma endregion




#pragma region WRITE TIME
/**
 * Write a time to a stream.
 * @param a the stream
 * @param x the time
 */
inline void writeTime(ostream& a, const time_t& x) {
  const int BUF = 64;
  char  buf[BUF];
  tm* t = localtime(&x);
  sprintf(buf, "%04d-%02d-%02d %02d:%02d:%02d",
    t->tm_year + 1900,
    t->tm_mon  + 1,
    t->tm_mday,
    t->tm_hour,
    t->tm_min,
    t->tm_sec
  );
  a << buf;
}


/**
 * Write a time point to a stream.
 * @param a the stream
 * @param x the time point
 */
inline void writeTimePoint(ostream& a, const time_point<system_clock>& x) {
  writeTime(a, system_clock::to_time_t(x));
}


/**
 * Write a time to a stream.
 * @param a the stream
 * @param x the time
 * @returns the stream
 */
inline ostream& operator<<(ostream& a, const time_t& x) {
  writeTime(a, x); return a;
}


/**
 * Write a time point to a stream.
 * @param a the stream
 * @param x the time point
 * @returns the stream
 */
inline ostream& operator<<(ostream& a, const time_point<system_clock>& x) {
  writeTimePoint(a, x); return a;
}
#pragma endregion




#pragma region PRINT
/**
 * Print an object.
 * @param x the object to print
 */
template <class T>
inline void print(const T& x) {
  cout << x;
}


/**
 * Print an object followed by a newline.
 * @param x the object to print
 */
template <class T>
inline void println(const T& x) {
  cout << x << "\n";
}


/**
 * Print a newline.
 */
inline void println() {
  cout << "\n";
}
#pragma endregion
