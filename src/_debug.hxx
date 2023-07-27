#pragma once
#include <chrono>
#include <ctime>
#include <cassert>
#include <cstdio>
#include <cstdlib>
#ifdef MPI
#include "_mpi.hxx"
#endif
#if !defined(NDEBUG) && defined(BUILD) && BUILD>=1
#include <unistd.h>
#include <signal.h>
#include <execinfo.h>
#endif

using std::chrono::system_clock;
using std::time_t;
using std::tm;




#pragma region BUILD MODES
#ifndef BUILD_RELEASE
/** Build has no debug information. */
#define BUILD_RELEASE 0
/** Build has only error information. */
#define BUILD_ERROR   1
/** Build has error and warning information. */
#define BUILD_WARNING 2
/** Build has error, warning and info information. */
#define BUILD_INFO    3
/** Build has error, warning, info and debug information. */
#define BUILD_DEBUG   4
/** Build has error, warning, info, debug and trace information. */
#define BUILD_TRACE   5
#endif
#pragma endregion




#pragma region PERFORM
#ifndef PEFORME
#if !defined(NDEBUG) && defined(BUILD) && BUILD>=BUILD_ERROR
/** Perform only if build mode is error or higher. */
#define PERFORME(...) __VA_ARGS__
#else
/** Perform only if build mode is error or higher. */
#define PERFORME(...)
#endif

#if !defined(NDEBUG) && defined(BUILD) && BUILD>=BUILD_WARNING
/** Perform only if build mode is warning or higher. */
#define PERFORMW(...) __VA_ARGS__
#else
/** Perform only if build mode is warning or higher. */
#define PERFORMW(...)
#endif

#if !defined(NDEBUG) && defined(BUILD) && BUILD>=BUILD_INFO
/** Perform only if build mode is info or higher. */
#define PERFORMI(...) __VA_ARGS__
#else
/** Perform only if build mode is info or higher. */
#define PERFORMI(...)
#endif

#if !defined(NDEBUG) && defined(BUILD) && BUILD>=BUILD_DEBUG
/** Perform only if build mode is debug or higher. */
#define PERFORMD(...) __VA_ARGS__
#else
/** Perform only if build mode is debug or higher. */
#define PERFORMD(...)
#endif

#if !defined(NDEBUG) && defined(BUILD) && BUILD>=BUILD_TRACE
/** Perform only if build mode is trace. */
#define PERFORMT(...) __VA_ARGS__
#else
/** Perform only if build mode is trace. */
#define PERFORMT(...)
#endif
#endif
#pragma endregion




#pragma region PRINT
#ifndef FPRINTFE
/** File print using format only if build mode is error or higher. */
#define FPRINTFE(...) PERFORME(fprintf(__VA_ARGS__))
/** File print using format only if build mode is warning or higher. */
#define FPRINTFW(...) PERFORMW(fprintf(__VA_ARGS__))
/** File print using format only if build mode is info or higher. */
#define FPRINTFI(...) PERFORMI(fprintf(__VA_ARGS__))
/** File print using format only if build mode is debug or higher. */
#define FPRINTFD(...) PERFORMD(fprintf(__VA_ARGS__))
/** File print using format only if build mode is trace. */
#define FPRINTFT(...) PERFORMT(fprintf(__VA_ARGS__))
#endif


#ifndef PRINTFE
/** Print using format only if build mode is error or higher. */
#define PRINTFE(...) PERFORME(printf(__VA_ARGS__))
/** Print using format only if build mode is warning or higher. */
#define PRINTFW(...) PERFORMW(printf(__VA_ARGS__))
/** Print using format only if build mode is info or higher. */
#define PRINTFI(...) PERFORMI(printf(__VA_ARGS__))
/** Print using format only if build mode is debug or higher. */
#define PRINTFD(...) PERFORMD(printf(__VA_ARGS__))
/** Print using format only if build mode is trace. */
#define PRINTFT(...) PERFORMT(printf(__VA_ARGS__))
#endif


#ifndef WRITEE
/** Write only if build mode is error or higher. */
#define WRITEE(...) PERFORME(write(__VA_ARGS__))
/** Write only if build mode is warning or higher. */
#define WRITEW(...) PERFORMW(write(__VA_ARGS__))
/** Write only if build mode is info or higher. */
#define WRITEI(...) PERFORMI(write(__VA_ARGS__))
/** Write only if build mode is debug or higher. */
#define WRITED(...) PERFORMD(write(__VA_ARGS__))
/** Write only if build mode is trace. */
#define WRITET(...) PERFORMT(write(__VA_ARGS__))
#endif


#ifndef PRINTE
/** Print only if build mode is error or higher. */
#define PRINTE(...) PERFORME(print(__VA_ARGS__))
/** Print only if build mode is warning or higher. */
#define PRINTW(...) PERFORMW(print(__VA_ARGS__))
/** Print only if build mode is info or higher. */
#define PRINTI(...) PERFORMI(print(__VA_ARGS__))
/** Print only if build mode is debug or higher. */
#define PRINTD(...) PERFORMD(print(__VA_ARGS__))
/** Print only if build mode is trace. */
#define PRINTT(...) PERFORMT(print(__VA_ARGS__))
#endif


#ifndef PRINTLNE
/** Print line only if build mode is error or higher. */
#define PRINTLNE(...) PERFORME(println(__VA_ARGS__))
/** Print line only if build mode is warning or higher. */
#define PRINTLNW(...) PERFORMW(println(__VA_ARGS__))
/** Print line only if build mode is info or higher. */
#define PRINTLNI(...) PERFORMI(println(__VA_ARGS__))
/** Print line only if build mode is debug or higher. */
#define PRINTLND(...) PERFORMD(println(__VA_ARGS__))
/** Print line only if build mode is trace. */
#define PRINTLNT(...) PERFORMT(println(__VA_ARGS__))
#endif
#pragma endregion




#pragma region LOG
#ifndef LOG
/**
 * Print log prefix.
 */
void logPrefix() {
  time_t s = system_clock::to_time_t(system_clock::now());
  tm    *t = localtime(&s);
#ifdef MPI
  printf("%04d-%02d-%02d %02d:%02d:%02d P%02d:"
#else
  printf("%04d-%02d-%02d %02d:%02d:%02d"
#endif
    , t->tm_year + 1900
    , t->tm_mon  + 1
    , t->tm_mday
    , t->tm_hour
    , t->tm_min
    , t->tm_sec
#ifdef MPI
    , mpi_comm_rank()
#endif
  );
}

/** Log using format. */
#define LOG(...) do { logPrefix(); printf(" " __VA_ARGS__); } while (0)
#endif

#ifndef LOGE
/** Log using format only if build mode is error or higher. */
#define LOGE(...) PERFORME(LOG(__VA_ARGS__))
/** Log using format only if build mode is warning or higher. */
#define LOGW(...) PERFORMW(LOG(__VA_ARGS__))
/** Log using format only if build mode is info or higher. */
#define LOGI(...) PERFORMI(LOG(__VA_ARGS__))
/** Log using format only if build mode is debug or higher. */
#define LOGD(...) PERFORMD(LOG(__VA_ARGS__))
/** Log using format only if build mode is trace. */
#define LOGT(...) PERFORMT(LOG(__VA_ARGS__))
#endif
#pragma endregion




#pragma region TRY
#ifndef TRY_MPIE
#ifdef  MPI
/** Try MPI expression only if build mode is error or higher. */
#define TRY_MPIE(exp) PERFORME(TRY_MPI(exp))
/** Try MPI expression only if build mode is warning or higher. */
#define TRY_MPIW(exp) PERFORMW(TRY_MPI(exp))
/** Try MPI expression only if build mode is info or higher. */
#define TRY_MPII(exp) PERFORMI(TRY_MPI(exp))
/** Try MPI expression only if build mode is debug or higher. */
#define TRY_MPID(exp) PERFORMD(TRY_MPI(exp))
/** Try MPI expression only if build mode is trace. */
#define TRY_MPIT(exp) PERFORMT(TRY_MPI(exp))
#endif
#endif


#ifndef TRY
#ifdef  MPI
/** Try expression. */
#define TRY(exp)  TRY_MPI(exp)
#endif
#endif


#ifndef TRYE
#ifdef  MPI
/** Try expression only if build mode is error or higher. */
#define TRYE(exp) TRY_MPIE(exp)
/** Try expression only if build mode is warning or higher. */
#define TRYW(exp) TRY_MPIW(exp)
/** Try expression only if build mode is info or higher. */
#define TRYI(exp) TRY_MPII(exp)
/** Try expression only if build mode is debug or higher. */
#define TRYD(exp) TRY_MPID(exp)
/** Try expression only if build mode is trace. */
#define TRYT(exp) TRY_MPIT(exp)
#endif
#endif
#pragma endregion




#pragma region ASSERT
#ifndef ASSERT
#if !defined(NDEBUG) && defined(BUILD) && BUILD>=BUILD_ERROR
#ifdef MPI
/** Assert expression. */
#define ASSERT(exp)           ASSERT_MPI(exp)
/** Assert expression with custom message. */
#define ASSERT_THAT(exp, msg) ASSERT_MPI((exp) && (msg))
#else
/** Assert expression. */
#define ASSERT(exp)           assert(exp)
/** Assert expression with custom message. */
#define ASSERT_THAT(exp, msg) assert((exp) && (msg))
#endif
#else
/** Assert expression. */
#define ASSERT(exp)
/** Assert expression with custom message. */
#define ASSERT_THAT(exp, msg)
#endif
#endif
#pragma endregion




#pragma region METHODS
#pragma region ON SIGNAL
/** Stack trace size for SIGSEGV signal handler. */
#define STACK_TRACE_SIZE 32

/**
 * Handle SIGSEGV signal.
 * @param sig signal number
 */
void on_sigsegv(int sig) {
#if !defined(NDEBUG) && defined(BUILD) && BUILD>=1
  void *entries[STACK_TRACE_SIZE];
  size_t n = backtrace(entries, STACK_TRACE_SIZE);
  fprintf(stderr, "ERROR: SIGNAL %d:\n", sig);
  backtrace_symbols_fd(entries, n, STDERR_FILENO);
  exit(1);
#endif
}
// - https://stackoverflow.com/a/77336/1413259


/**
 * Install SIGSEGV signal handler, if build mode is error or higher.
 */
void install_sigsegv() {
#if !defined(NDEBUG) && defined(BUILD) && BUILD>=1
  signal(SIGSEGV, on_sigsegv);
#endif
}
// - https://stackoverflow.com/a/77336/1413259
#pragma endregion
#pragma endregion
