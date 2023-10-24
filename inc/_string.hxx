#pragma once
#include <string>
#include <cstdint>
#include "_debug.hxx"

using std::string;




#pragma region METHODS
/**
 * Count the number of lines in a string.
 * @param x string
 * @returns number of lines
 */
inline size_t countLines(const char* x) {
  ASSERT(x);
  size_t a = 1;
  for (; *x; x++) {
    if (*x == '\r' || *x == '\n') ++a;
    else if (*x == '\r' && *(x+1) == '\n') ++x;
  }
  return a;
}


/**
 * Count the number of lines in a string.
 * @param x string
 * @returns number of lines
 */
inline size_t countLines(const string& x) {
  return countLines(x.c_str());
}
#pragma endregion
