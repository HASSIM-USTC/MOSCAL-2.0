#pragma once

#include "deom.hpp"

inline int INDEX2(const int i, const int j, const int ld) { return i * ld + j; }
inline int INDEX3(const int i, const int j, const int k, const int ld) {
  return i * ld * ld + j * ld + k;
}
inline int INDEX2(const int i, const int j) { return i * NSYS + j; }
inline int INDEX3(const int i, const int j, const int k) {
  return i * NSYS * NSYS + j * NSYS + k;
}