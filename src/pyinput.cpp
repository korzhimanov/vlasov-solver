/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file pyinput.cpp
 * \brief The source file which implements the methods of the pyinput class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#include "include/pyinput.h"
#include <iostream>

#include "include/errors.h"

pyinput::pyinput(int &err) {
  Py_Initialize();
  if (!Py_IsInitialized()) {
    PyErr_Print();
    std::cerr << "Terminated!" << std::endl;
    err = PYTHON_INITIALIZATION_FAILED;
    return;
  }
  main_module = PyImport_AddModule("__main__");
  main_dict = PyModule_GetDict(main_module);
  PyRun_SimpleString("from math import *");
  PyRun_SimpleString("def sqr(x):\n\treturn x*x");
}

pyinput::~pyinput() {
  if (Py_IsInitialized()) Py_Finalize();
}

/**
 * \todo Add error output
 */
void pyinput::ReadFile(std::string filename, int &err) {
  FILE *f = fopen(filename.c_str(), "r");
  if (f == NULL) {
    std::cerr << "Cannot open init-file " << filename << "." << std::endl;
    err = CANNOT_OPEN_INIT_FILE;
    return;
  } else {
    err = PyRun_SimpleFile(f, filename.c_str());
    fclose(f);
    if (err != 0) {
      std::cerr << "There was an error in parsing init-file. Be sure it's a "
                   "correct python script." << std::endl;
    }
  }
}

void pyinput::Set(std::string name, int &var, int &err) const {
  PyObject *io = GetVarObject(name, err);
  if (err == 0) {
    if (!PyInt_Check(io)) {
      std::cerr << name << " must be integer." << std::endl;
      err = VAR_MUST_INTEGER;
    }
    var = PyInt_AsLong(io);
  }
}

void pyinput::Set(std::string name, double &var, int &err) const {
  PyObject *io = GetVarObject(name, err);
  if (err == 0) {
    if (!PyFloat_Check(io)) {
      std::cerr << name << " must be float." << std::endl;
      err = VAR_MUST_FLOAT;
    }
    var = PyFloat_AsDouble(io);
  }
}

void pyinput::Set(std::string name, bool &var, int &err) const {
  PyObject *io = GetVarObject(name, err);
  if (err == 0) {
    if (!PyBool_Check(io)) {
      std::cerr << name << " must be boolean." << std::endl;
      err = VAR_MUST_BOOLEAN;
    }
    var = PyInt_AsLong(io);
  }
}

void pyinput::Set(std::string name, std::string &var, int &err) const {
  PyObject *io = GetVarObject(name, err);
  if (err == 0) {
    if (!PyString_Check(io)) {
      std::cerr << name << " must be string." << std::endl;
      err = VAR_MUST_STRING;
    }
    var = std::string(PyString_AsString(io));
  }
}

void pyinput::SetNotNegative(std::string name, int &var, int &err) const {
  Set(name, var, err);
  if (err) return;
  if (var < 0) {
    std::cerr << name + " mustnot be negative." << std::endl;
    err = VAR_MUST_NOTNEGATIVE;
  }
}

void pyinput::SetPositive(std::string name, int &var, int &err) const {
  Set(name, var, err);
  if (err) return;
  if (var <= 0) {
    std::cerr << name + " must be positive." << std::endl;
    err = VAR_MUST_POSITIVE;
  }
}

void pyinput::SetPositive(std::string name, double &var, int &err) const {
  Set(name, var, err);
  if (err) return;
  if (var <= 0.) {
    std::cerr << name + " must be positive." << std::endl;
    err = VAR_MUST_POSITIVE;
  }
}

pFunc pyinput::GetFunc(std::string name, int &err) const {
  PyObject *io = GetVarObject(name, err);
  if (err != 0) {
    return pFunc();
  }
  if (err == 0) {
    if (!PyFunction_Check(io)) {
      std::cerr << name << " must be function." << std::endl;
      err = VAR_MUST_FUNCTION;
      return pFunc();
    }
    return pFunc(io);
  }
}

PyObject *pyinput::GetVarObject(std::string &name, int &err) const {
  PyObject *io = PyDict_GetItemString(main_dict, name.c_str());
  if (io == NULL) {
    std::cerr << "Parameter " << name << " not found in the init-file."
              << std::endl;
    err = VAR_NOT_FOUND;
  }
  return io;
}
