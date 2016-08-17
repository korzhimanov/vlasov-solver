/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file mymath.h
 * \brief The header file which defines MyMath class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#ifndef MYMATH_H
#define MYMATH_H

#include <cstring>
#include <cmath>
#include <gsl/gsl_math.h>

/**
 * /class MyMath
 * The class implements some useful mathematical fuctions
*/
class MyMath
{
public:
    MyMath();
    virtual ~MyMath();
    void zeros(double *a, const int number);
    void zeros(double *a, const int number, const int step);

    double sum(double *a, const int number, const int step = 1);

    template<typename Type> Type max(Type *a, const int number, const int step = 1);
    template<typename Type> Type min(Type *a, const int number, const int step = 1);
    template<typename Type> Type sqr(Type x);
};

inline void MyMath::zeros(double *a, const int number)
{
    memset(a, 0, sizeof(*a)*number);
}

template<typename Type>
Type MyMath::max(Type *a, const int number, const int step)
{
    Type tmp = *a, *t = a+1, *until = a + number*step;
    for (; t < until; t += step) tmp = (tmp < *t) ? *t : tmp;
    return tmp;
}

template<typename Type>
Type MyMath::min(Type *a, const int number, const int step)
{
    Type tmp = *a, *t = a+1, *until = a + number*step;
    for (; t < until; t += step) tmp = (tmp > *t) ? *t : tmp;
    return tmp;
}

template<typename Type>
inline Type MyMath::sqr(Type x)
{
    return x*x;
}

#endif // MYMATH_H
