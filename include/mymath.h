/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file mymath.h
 * \brief The header file which defines some mathematical functions
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#ifndef INCLUDE_MYMATH_H_
#define INCLUDE_MYMATH_H_

#include <cmath>
#include <cstring>

namespace mymath {

void zeros(double *a, const int number);
void zeros(double *a, const int number, const int step);

double sum(double *a, const int number, const int step = 1);

template <typename Type>
Type sqr(Type x);

double const onesixth = 1. / 6.;

}  // namespace mymath

inline void mymath::zeros(double *a, const int number) {
  memset(a, 0, sizeof(*a) * number);
}

template <typename Type>
inline Type mymath::sqr(Type x) {
  return x * x;
}

#endif  // INCLUDE_MYMATH_H_
