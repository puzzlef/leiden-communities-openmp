#pragma once
#include <iterator>
#include <cstdint>

using std::iterator_traits;




#pragma region CLASSES
/**
 * A bounded (circular) deque view capable of storing N elements.
 * @tparam I iterator type of the buffer used to store the elements
 */
template <class I>
class DequeView {
  #pragma region TYPES
  protected:
  using T = typename iterator_traits<I>::value_type;

  public:
  /** Value type of the deque. */
  using value_type = T;
  #pragma endregion


  #pragma region DATA
  protected:
  /** Begin iterator of the buffer. */
  const I xb;
  /** End iterator of the buffer. */
  const I xe;
  /** Front iterator of the deque. */
  I ib;
  /** Rear iterator of the deque. */
  I ie;
  /** Number of elements in the deque. */
  size_t n;
  #pragma endregion


  #pragma region METHODS
  #pragma region SIZE
  public:
  /**
   * Get the size of the deque.
   * @returns the size
   */
  inline size_t size() const {
    return n;
  }

  /**
   * Check if the deque is empty.
   * @returns true if empty
   */
  inline bool empty() const {
    return n==0;
  }
  #pragma endregion


  #pragma region READ
  public:
  /**
   * Get the element at the rear of the deque.
   * @returns the element
   */
  inline auto back() const {
    return ie==xb? *(xe-1) : *(ie-1);
  }

  /**
   * Get the element at the front of the deque.
   * @returns the element
   */
  inline auto front() const {
    return *ib;
  }
  #pragma endregion


  #pragma region WRITE
  public:
  /**
   * Push value to the rear of the deque.
   * @param v value to push
   */
  inline void push_back(const T& v) {
    if (ib!=ie || n==0) ++n;
    *(ie++) = v;
    if (ie==xe) ie = xb;
  }

  /**
   * Push value to the front of the deque.
   * @param v value to push
   */
  inline void push_front(const T& v) {
    if (ib==xb) ib = xe;
    *(--ib) = v;
    if (ib!=ie || n==0) ++n;
  }

  /**
   * Pop value from the rear of the deque.
   * @returns the value
   */
  inline auto pop_back() {
    if (ie==xb) ie = xe;
    if (n>0) --n;
    return *(--ie);
  }

  /**
   * Pop value from the front of the deque.
   * @returns the value
   */
  inline auto pop_front() {
    if (n>0) --n;
    auto v = *(ib++);
    if (ib==xe) ib = xb;
    return v;
  }
  #pragma endregion
  #pragma endregion


  #pragma region CONSTRUCTORS
  public:
  /**
   * Construct a Deque View from a given buffer.
   * @param xb begin iterator of the buffer
   * @param xe end iterator of the buffer
   */
  DequeView(I xb, I xe) :
  xb(xb), xe(xe), ib(xb), ie(xb), n(0) {}
  #pragma endregion
};


/**
 * Construct a deque view from a given buffer.
 * @param xb begin iterator of the buffer
 * @param xe end iterator of the buffer
 * @returns the Deque View
 */
template <class I>
inline auto deque_view(I xb, I xe) {
  return DequeView<I>(xb, xe);
}




/**
 * A bounded (circular) deque view capable of storing N-1 elements.
 * @tparam I iterator type of the buffer used to store the elements
 */
template <class I>
class UnsizedDequeView {
  #pragma region TYPES
  protected:
  using T = typename iterator_traits<I>::value_type;

  public:
  /** Value type of the deque. */
  using value_type = T;
  #pragma endregion


  #pragma region DATA
  protected:
  /** Begin iterator of the buffer. */
  const I xb;
  /** End iterator of the buffer. */
  const I xe;
  /** Front iterator of the deque. */
  I ib;
  /** Rear iterator of the deque. */
  I ie;
  #pragma endregion


  #pragma region METHODS
  #pragma region SIZE
  public:
  /**
   * Check if the deque is empty.
   * @returns true if empty
   */
  inline bool empty() const {
    return ib==ie;
  }
  #pragma endregion


  #pragma region READ
  public:
  /**
   * Get the element at the rear of the deque.
   * @returns the element
   */
  inline auto back()  const {
    return ie==xb? *(xe-1) : *(ie-1);
  }

  /**
   * Get the element at the front of the deque.
   * @returns the element
   */
  inline auto front() const {
    return *ib;
  }
  #pragma endregion


  #pragma region WRITE
  public:
  /**
   * Push value to the rear of the deque.
   * @param v value to push
   */
  inline void push_back(const T& v) {
    *(ie++) = v;
    if (ie==xe) ie = xb;
  }

  /**
   * Push value to the front of the deque.
   * @param v value to push
   */
  inline void push_front(const T& v) {
    if (ib==xb) ib = xe;
    *(--ib) = v;
  }

  /**
   * Pop value from the rear of the deque.
   * @returns the value
   */
  inline auto pop_back() {
    if (ie==xb) ie = xe;
    return *(--ie);
  }

  /**
   * Pop value from the front of the deque.
   * @returns the value
   */
  inline auto pop_front() {
    auto v = *(ib++);
    if (ib==xe) ib = xb;
    return v;
  }
  #pragma endregion
  #pragma endregion


  #pragma region CONSTRUCTORS
  public:
  /**
   * Construct an unsized deque view from a given buffer.
   * @param xb begin iterator of the buffer
   * @param xe end iterator of the buffer
   */
  UnsizedDequeView(I xb, I xe) :
  xb(xb), xe(xe), ib(xb), ie(xb) {}
  #pragma endregion
};


/**
 * Construct an unsized deque view from a given buffer.
 * @param xb begin iterator of the buffer
 * @param xe end iterator of the buffer
 * @returns the deque view
 */
template <class I>
inline auto unsized_deque_view(I xb, I xe) {
  return UnsizedDequeView<I>(xb, xe);
}
#pragma endregion
