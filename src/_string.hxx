#pragma once
#include <string>
#include "_debug.hxx"

using std::string;




// COUNT LINES
// -----------
// For counting temporal edges.

size_t countLines(const char* x) {
  ASSERT(x);
  size_t a = 1;
  for (; *x; x++) {
    if (*x == '\r' || *x == '\n') ++a;
    else if (*x == '\r' && *(x+1) == '\n') ++x;
  }
  return a;
}
inline size_t countLines(const string& x) {
  return countLines(x.c_str());
}
