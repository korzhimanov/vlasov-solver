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

#ifndef INCLUDE_PFUNC_H_
#define INCLUDE_PFUNC_H_

#include <Python.h>

/**
 * /class pFunc
 * The class for using functions defined in python scripts
*/
class pFunc {
 public:
  explicit pFunc(PyObject *);
  pFunc(const pFunc &);
  virtual ~pFunc();
  double call(double);
  double call(double, double);
  double call(double, double, double);

 public:
  PyObject *function;
};

#endif  // INCLUDE_PFUNC_H_
