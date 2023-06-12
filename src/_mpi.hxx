#pragma once
#include <vector>
#include <mpi.h>




// TRY
// ---
// Log error if MPI function call fails.

#ifndef TRY_MPI
void tryFailedMpi(int err, const char* exp, const char* func, int line, const char* file) {
  char buf[MPI_MAX_ERROR_STRING]; int BUF;
  MPI_Error_string(err, buf, &BUF);
  fprintf(stderr,
    "ERROR: %s\n"
    "  in expression %s\n"
    "  at %s:%d in %s\n",
    buf, exp, func, line, file);
  MPI_Abort(MPI_COMM_WORLD, err);
}
#define TRY_MPI(exp)  do { int err = exp; if (err != MPI_SUCCESS) tryFailedMpi(err, #exp, __func__, __LINE__, __FILE__); } while (0)
#endif




// ASSERT
// ------
// Assert that expression is true.

#ifndef ASSERT_MPI
void assertFailedMpi(const char* exp, const char* func, int line, const char* file) {
  fprintf(stderr,
    "ERROR: Assertion failed\n"
    "  in expression %s\n"
    "  at %s:%d in %s\n",
    exp, func, line, file);
  MPI_Abort(MPI_COMM_WORLD, 1);
}
#define ASSERT_MPI(exp)  do { if (!(exp)) assertFailedMpi(#exp, __func__, __LINE__, __FILE__); } while (0)
#endif




// BASIC
// -----

inline int mpi_comm_size(MPI_Comm comm=MPI_COMM_WORLD) {
  int size; MPI_Comm_size(comm, &size);
  return size;
}

inline int mpi_comm_rank(MPI_Comm comm=MPI_COMM_WORLD) {
  int rank; MPI_Comm_rank(comm, &rank);
  return rank;
}
