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




// BUILD
// -----
// Build modes.

#ifndef BUILD_RELEASE
#define BUILD_RELEASE 0
#define BUILD_ERROR   1
#define BUILD_WARNING 2
#define BUILD_INFO    3
#define BUILD_DEBUG   4
#define BUILD_TRACE   5
#endif




// PERFORM
// -------

#ifndef PEFORME
#if !defined(NDEBUG) && defined(BUILD) && BUILD>=BUILD_ERROR
#define PERFORME(...) __VA_ARGS__
#else
#define PERFORME(...)
#endif
#if !defined(NDEBUG) && defined(BUILD) && BUILD>=BUILD_WARNING
#define PERFORMW(...) __VA_ARGS__
#else
#define PERFORMW(...)
#endif
#if !defined(NDEBUG) && defined(BUILD) && BUILD>=BUILD_INFO
#define PERFORMI(...) __VA_ARGS__
#else
#define PERFORMI(...)
#endif
#if !defined(NDEBUG) && defined(BUILD) && BUILD>=BUILD_DEBUG
#define PERFORMD(...) __VA_ARGS__
#else
#define PERFORMD(...)
#endif
#if !defined(NDEBUG) && defined(BUILD) && BUILD>=BUILD_TRACE
#define PERFORMT(...) __VA_ARGS__
#else
#define PERFORMT(...)
#endif
#endif




// PRINT
// -----

#ifndef FPRINTFE
#define FPRINTFE(...) PERFORME(fprintf(__VA_ARGS__))
#define FPRINTFW(...) PERFORMW(fprintf(__VA_ARGS__))
#define FPRINTFI(...) PERFORMI(fprintf(__VA_ARGS__))
#define FPRINTFD(...) PERFORMD(fprintf(__VA_ARGS__))
#define FPRINTFT(...) PERFORMT(fprintf(__VA_ARGS__))
#endif

#ifndef PRINTFE
#define PRINTFE(...) PERFORME(printf(__VA_ARGS__))
#define PRINTFW(...) PERFORMW(printf(__VA_ARGS__))
#define PRINTFI(...) PERFORMI(printf(__VA_ARGS__))
#define PRINTFD(...) PERFORMD(printf(__VA_ARGS__))
#define PRINTFT(...) PERFORMT(printf(__VA_ARGS__))
#endif

#ifndef WRITEE
#define WRITEE(...) PERFORME(write(__VA_ARGS__))
#define WRITEW(...) PERFORMW(write(__VA_ARGS__))
#define WRITEI(...) PERFORMI(write(__VA_ARGS__))
#define WRITED(...) PERFORMD(write(__VA_ARGS__))
#define WRITET(...) PERFORMT(write(__VA_ARGS__))
#endif

#ifndef PRINTE
#define PRINTE(...) PERFORME(print(__VA_ARGS__))
#define PRINTW(...) PERFORMW(print(__VA_ARGS__))
#define PRINTI(...) PERFORMI(print(__VA_ARGS__))
#define PRINTD(...) PERFORMD(print(__VA_ARGS__))
#define PRINTT(...) PERFORMT(print(__VA_ARGS__))
#endif

#ifndef PRINTLNE
#define PRINTLNE(...) PERFORME(println(__VA_ARGS__))
#define PRINTLNW(...) PERFORMW(println(__VA_ARGS__))
#define PRINTLNI(...) PERFORMI(println(__VA_ARGS__))
#define PRINTLND(...) PERFORMD(println(__VA_ARGS__))
#define PRINTLNT(...) PERFORMT(println(__VA_ARGS__))
#endif




// LOG
// ---

#ifndef LOG
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
#define LOG(...) do { logPrefix(); printf(" " __VA_ARGS__); } while (0)
#endif

#ifndef LOGE
#define LOGE(...) PERFORME(LOG(__VA_ARGS__))
#define LOGW(...) PERFORMW(LOG(__VA_ARGS__))
#define LOGI(...) PERFORMI(LOG(__VA_ARGS__))
#define LOGD(...) PERFORMD(LOG(__VA_ARGS__))
#define LOGT(...) PERFORMT(LOG(__VA_ARGS__))
#endif




// TRY
// ---

#ifndef TRY_MPIE
#ifdef  MPI
#define TRY_MPIE(exp) PERFORME(TRY_MPI(exp))
#define TRY_MPIW(exp) PERFORMW(TRY_MPI(exp))
#define TRY_MPII(exp) PERFORMI(TRY_MPI(exp))
#define TRY_MPID(exp) PERFORMD(TRY_MPI(exp))
#define TRY_MPIT(exp) PERFORMT(TRY_MPI(exp))
#endif
#endif


#ifndef TRY
#ifdef  MPI
#define TRY(exp)  TRY_MPI(exp)
#endif
#endif

#ifndef TRYE
#ifdef  MPI
#define TRYE(exp) TRY_MPIE(exp)
#define TRYW(exp) TRY_MPIW(exp)
#define TRYI(exp) TRY_MPII(exp)
#define TRYD(exp) TRY_MPID(exp)
#define TRYT(exp) TRY_MPIT(exp)
#endif
#endif




// ASSERT
// ------

#ifndef ASSERT
#if !defined(NDEBUG) && defined(BUILD) && BUILD>=BUILD_ERROR
#ifdef MPI
#define ASSERT(exp)           ASSERT_MPI(exp)
#define ASSERT_THAT(exp, msg) ASSERT_MPI((exp) && (msg))
#else
#define ASSERT(exp)           assert(exp)
#define ASSERT_THAT(exp, msg) assert((exp) && (msg))
#endif
#else
#define ASSERT(exp)
#define ASSERT_THAT(exp, msg)
#endif
#endif




// ON SIGNAL
// ---------

#define STACK_TRACE_SIZE 32

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

void install_sigsegv() {
#if !defined(NDEBUG) && defined(BUILD) && BUILD>=1
  signal(SIGSEGV, on_sigsegv);
#endif
}
// - https://stackoverflow.com/a/77336/1413259
