#pragma once
#include <cstdint>




#pragma region CLASSES
/**
 * A 32-bit xorshift RNG.
 */
class xorshift32_engine {
  #pragma region DATA
  private:
  /** State of the RNG. */
  uint32_t state;
  #pragma endregion


  #pragma region METHODS
  public:
  /**
   * Generate a random number.
   */
  uint32_t operator()() {
    uint32_t x = state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    return state = x;
  }
  #pragma endregion


  #pragma region CONSTRUCTORS
  /**
   * Construct an RNG with a random seed.
   * @param state initial state
   */
  xorshift32_engine(uint32_t state)
  : state(state) {}
  #pragma endregion
};
// - https://stackoverflow.com/a/71523041/1413259
// - https://www.jstatsoft.org/article/download/v008i14/916
#pragma endregion
