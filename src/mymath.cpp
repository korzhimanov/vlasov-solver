/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file mymath.cpp
 * \brief The source file which implements some mathematical functions
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#include "mymath.h"

namespace mymath {

void zeros(double *a, const int number, const int step)
{
    double *t = a, *until = a + number*step;
    for (; t < until; t += step) *t = 0.;
}

double sum(double *a, const int number, const int step)
{
    double tmp = *a, *t = a+1, *until = a + number*step;
    for (; t < until; t += step) tmp += *t;
    return tmp;
}

} // namespace mymath
