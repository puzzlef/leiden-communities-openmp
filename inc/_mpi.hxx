#pragma once
#include <chrono>
#include <ctime>
#include <cstdio>
#include <cstdlib>
#include <mpi.h>
#include "_debug.hxx"

using std::chrono::system_clock;
using std::time_t;
using std::tm;
using std::localtime;
using std::fprintf;
using std::printf;




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

/**
 * Try to execute an MPI function call only if build mode is error or higher.
 * @param exp expression to execute
 **/
#define TRY_MPIE(exp)  PERFORME(TRY_MPI(exp))

/**
 * Try to execute an MPI function call only if build mode is warning or higher.
 * @param exp expression to execute
 **/
#define TRY_MPIW(exp)  PERFORMW(TRY_MPI(exp))

/**
 * Try to execute an MPI function call only if build mode is info or higher.
 * @param exp expression to execute
 **/
#define TRY_MPII(exp)  PERFORMI(TRY_MPI(exp))

/**
 * Try to execute an MPI function call only if build mode is debug or higher.
 * @param exp expression to execute
 **/
#define TRY_MPID(exp)  PERFORMD(TRY_MPI(exp))

/**
 * Try to execute an MPI function call only if build mode is trace.
 * @param exp expression to execute
 **/
#define TRY_MPIT(exp)  PERFORMT(TRY_MPI(exp))
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




#pragma region LOG
#ifndef LOG_MPI
/**
 * Print log prefix.
 */
void logPrefixMpi() {
  time_t s = system_clock::to_time_t(system_clock::now());
  tm    *t = localtime(&s);
  printf("%04d-%02d-%02d %02d:%02d:%02d P%02d:"
    , t->tm_year + 1900
    , t->tm_mon  + 1
    , t->tm_mday
    , t->tm_hour
    , t->tm_min
    , t->tm_sec
    , mpi_comm_rank()
  );
}

/** Log using format. */
#define LOG_MPI(...) do { logPrefixMpi(); printf(" " __VA_ARGS__); } while (0)
#endif
#pragma endregion
