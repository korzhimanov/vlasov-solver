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

#include "pyinput.h"
#include <iostream>

pyinput::pyinput()
{
    Py_Initialize();
    main_module = PyImport_AddModule("__main__");
    main_dict = PyModule_GetDict(main_module);
    PyRun_SimpleString("from math import *");
    PyRun_SimpleString("def sqr(x):\n\treturn x*x");
}

pyinput::~pyinput()
{
    if (f != NULL) fclose(f);
    Py_Finalize();
}

/**
 * \todo Add error output
 */
void pyinput::ReadFile(std::string filename)
{
    f = fopen(filename.c_str(), "r");
    PyRun_SimpleFile(f, filename.c_str());
}

int pyinput::GetInt(std::string name)
{
    return PyInt_AsLong(PyDict_GetItemString(main_dict, name.c_str()));
}

double pyinput::GetDouble(std::string name)
{
    return PyFloat_AsDouble(PyDict_GetItemString(main_dict, name.c_str()));
}

std::string pyinput::GetString(std::string name)
{
    return std::string(PyString_AsString(PyDict_GetItemString(main_dict, name.c_str())));
}

pFunc pyinput::GetFunc(std::string name)
{
    pFunc func(PyDict_GetItemString(main_dict, name.c_str()));
    return func;
}

bool pyinput::Set(std::string name, int *var)
{
    *var = GetInt(name);
    return true;
}

bool pyinput::Set(std::string name, double *var)
{
    *var = GetDouble(name);
    return true;
}

bool pyinput::SetNotNegative(std::string name, int *var)
{
    *var = GetInt(name);
    if ( *var >= 0 ) return true;
    else {std::cout << name + " mustnot be negative" << std::endl; return false;}
}

bool pyinput::SetPositive(std::string name, int *var)
{
    *var = GetInt(name);
    if ( *var > 0 ) return true;
    else {std::cout << name + " must be positive" << std::endl; return false;}
}

bool pyinput::SetPositive(std::string name, double *var)
{
    *var = GetDouble(name);
    if ( *var > 0. ) return true;
    else {std::cout << name + " must be positive" << std::endl; return false;}
}
