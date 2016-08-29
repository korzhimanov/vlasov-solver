/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file pyfunc.h
 * \brief The header file which defines PyFunc class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#ifndef INCLUDE_PFUNC_H_
#define INCLUDE_PFUNC_H_

#include <Python.h>

/**
 * /class PyFunc
 * The class for using functions defined in python scripts
*/
class PyFunc {
 public:
  PyFunc();
  explicit PyFunc(PyObject *);
  PyFunc(const PyFunc &);
  virtual ~PyFunc();
  double call(double);
  double call(double, double);
  double call(double, double, double);

 public:
  PyObject *function;
};

#endif  // INCLUDE_PFUNC_H_
