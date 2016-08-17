/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file pfunc.h
 * \brief The header file which defines pFunc class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#ifndef PFUNC_H
#define PFUNC_H

#include <Python.h>

/**
 * /class pFunc
 * The class for using functions defined in python scripts
*/
class pFunc
{
    public:
        pFunc(PyObject*);
        pFunc(const pFunc&);
        virtual ~pFunc();
        double call(double);
        double call(double, double);
        double call(double, double, double);

    public:
        PyObject *function;
};

#endif // PFUNC_H
