#pragma once
#include <type_traits>
#include <istream>
#include <ostream>
#include <cstdint>
#include "_debug.hxx"

using std::make_signed_t;
using std::istream;
using std::ostream;




#pragma region TYPES
/** Signed size type. */
using ssize_t = make_signed_t<size_t>;
#pragma endregion




#pragma region CLASSES
#ifndef NONE
/**
 * Zero size type.
 */
struct None {
  #pragma region METHODS
  #pragma region COMPARISION OPERATORS
  friend inline bool operator==(None l, None r)     noexcept { return true; }
  template <class T>
  friend inline bool operator==(None l, const T& r) noexcept { return false; }
  template <class T>
  friend inline bool operator==(const T& l, None r) noexcept { return false; }
  #pragma endregion


  #pragma region STREAM OPERATORS
  friend inline istream& operator>>(istream& a, None& x) noexcept { return a; }
  friend inline ostream& operator<<(ostream& a, None x)  noexcept { return a; }
  #pragma endregion
  #pragma endregion


  #pragma region CONSTRUCTORS
  /**
   * Construct a zero size object.
   */
  explicit None() {}

  /**
   * Construct a zero size object.
   * @param _ any value (ignored)
   */
  template <class T>
  explicit None(T _) {}
  #pragma endregion
};
#define NONE None
#endif




#ifndef SHOW_TYPE
/**
 * Show type of template parameter at compile time.
 * @tparam T type to show
 */
template <class T>
struct ShowType;

/**
 * Show type of template parameter at compile time.
 * @param T type to show
 */
#define SHOW_TYPE(T) ShowType<T>{}

/**
 * Show type of template parameter at compile time, only if build mode is error or higher.
 * @param T type to show
 */
#define SHOW_TYPEE(T)  PERFORME(SHOW_TYPE(T))

/**
 * Show type of template parameter at compile time, only if build mode is warning or higher.
 * @param T type to show
 */
#define SHOW_TYPEW(T)  PERFORMW(SHOW_TYPE(T))

/**
 * Show type of template parameter at compile time, only if build mode is info or higher.
 * @param T type to show
 */
#define SHOW_TYPEI(T)  PERFORMI(SHOW_TYPE(T))

/**
 * Show type of template parameter at compile time, only if build mode is debug or higher.
 * @param T type to show
 */
#define SHOW_TYPED(T)  PERFORMD(SHOW_TYPE(T))

/**
 * Show type of template parameter at compile time, only if build mode is trace or higher.
 * @param T type to show
 */
#define SHOW_TYPET(T)  PERFORMT(SHOW_TYPE(T))
#endif
#pragma endregion
// - https://stackoverflow.com/a/39637476/1413259




#pragma region TEMPLATE TYPES
// For template classes in template :).
#ifndef tclass0
/** Template class with no arguments. */
#define tclass0  template <> class
/** Template class with one argument. */
#define tclass1  template <class> class
/** Template class with two arguments. */
#define tclass2  template <class, class> class
/** Template class with three arguments. */
#define tclass3  template <class, class, class> class
/** Template class with zero or more arguments. */
#define tclass0s  template <class...> class
/** Template class with one or more arguments. */
#define tclass1s  template <class, class...> class
/** Template class with two or more arguments. */
#define tclass2s  template <class, class, class...> class
/** Template class with three or more arguments. */
#define tclass3s  template <class, class, class, class...> class
#endif
#pragma endregion
