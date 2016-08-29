/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file pyinput.h
 * \brief The header file which defines PyInput class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#ifndef INCLUDE_PYINPUT_H_
#define INCLUDE_PYINPUT_H_

#include <string>

#include "include/pyfunc.h"

/**
 * /class PyFunc
 * The class for reading python script with input parameters
*/
class PyInput {
 public:
  PyInput(int &err);
  virtual ~PyInput();

  void ReadFile(std::string, int &err);

  void Set(std::string name, int &var, int &err) const;
  void Set(std::string name, double &var, int &err) const;
  void Set(std::string name, bool &var, int &err) const;
  void Set(std::string name, std::string &var, int &err) const;
  void SetNotNegative(std::string name, int &var, int &err) const;
  void SetPositive(std::string name, int &var, int &err) const;
  void SetPositive(std::string name, double &var, int &err) const;

  PyFunc GetFunc(std::string, int &err) const;

 private:
  PyObject *GetVarObject(std::string &, int &err) const;

 private:
  PyObject *main_module, *main_dict;
};

#endif  // INCLUDE_PYINPUT_H_
