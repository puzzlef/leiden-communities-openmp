#pragma once
#include <mpi.h>




#pragma region METHODS
#pragma region ERROR
#ifndef TRY_MPI
/**
 * Log error on MPI function call failure.
 * @param err error code
 * @param exp expression string
 * @param func current function name
 * @param line current line number
 * @param file current file name
 */
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

/**
 * Try to execute an MPI function call.
 * @param exp expression to execute
 */
#define TRY_MPI(exp)  do { int err = exp; if (err != MPI_SUCCESS) tryFailedMpi(err, #exp, __func__, __LINE__, __FILE__); } while (0)
#endif




#ifndef ASSERT_MPI
/**
 * Log error on assertion failure.
 * @param exp expression string
 * @param func current function name
 * @param line current line number
 * @param file current file name
 */
void assertFailedMpi(const char* exp, const char* func, int line, const char* file) {
  fprintf(stderr,
    "ERROR: Assertion failed\n"
    "  in expression %s\n"
    "  at %s:%d in %s\n",
    exp, func, line, file);
  MPI_Abort(MPI_COMM_WORLD, 1);
}

/**
 * Assert that expression is true.
 * @param exp expression that should be true
 */
#define ASSERT_MPI(exp)  do { if (!(exp)) assertFailedMpi(#exp, __func__, __LINE__, __FILE__); } while (0)
#endif
#pragma endregion




#pragma region BASIC
/**
 * Get the size of the communicator.
 * @param comm communicator
 * @returns the size
 */
inline int mpi_comm_size(MPI_Comm comm=MPI_COMM_WORLD) {
  int size; MPI_Comm_size(comm, &size);
  return size;
}


/**
 * Get the rank of the current process.
 * @param comm communicator
 * @returns the rank
 */
inline int mpi_comm_rank(MPI_Comm comm=MPI_COMM_WORLD) {
  int rank; MPI_Comm_rank(comm, &rank);
  return rank;
}
#pragma endregion
#pragma endregion
