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

#ifndef MYMATH_H
#define MYMATH_H

#include <cstring>
#include <cmath>
#include <gsl/gsl_math.h>

namespace mymath {

void zeros(double *a, const int number);
void zeros(double *a, const int number, const int step);

double sum(double *a, const int number, const int step = 1);

template<typename Type> Type max(Type *a, const int number, const int step = 1);
template<typename Type> Type min(Type *a, const int number, const int step = 1);
template<typename Type> Type sqr(Type x);

} // namespace mymath

inline void mymath::zeros(double *a, const int number)
{
    memset(a, 0, sizeof(*a)*number);
}

template<typename Type>
Type mymath::max(Type *a, const int number, const int step)
{
    Type tmp = *a, *t = a+1, *until = a + number*step;
    for (; t < until; t += step) tmp = (tmp < *t) ? *t : tmp;
    return tmp;
}

template<typename Type>
Type mymath::min(Type *a, const int number, const int step)
{
    Type tmp = *a, *t = a+1, *until = a + number*step;
    for (; t < until; t += step) tmp = (tmp > *t) ? *t : tmp;
    return tmp;
}

template<typename Type>
inline Type mymath::sqr(Type x)
{
    return x*x;
}

#endif // MYMATH_H
