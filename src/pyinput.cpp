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

using namespace std;

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

void pyinput::ReadFile(string filename)
{
    f = fopen(filename.c_str(), "r");
    PyRun_SimpleFile(f, filename.c_str());
}

int pyinput::GetInt(string name)
{
    return PyInt_AsLong(PyDict_GetItemString(main_dict, name.c_str()));
}

double pyinput::GetDouble(string name)
{
    return PyFloat_AsDouble(PyDict_GetItemString(main_dict, name.c_str()));
}

string pyinput::GetString(string name)
{
    return string(PyString_AsString(PyDict_GetItemString(main_dict, name.c_str())));
}

pFunc pyinput::GetFunc(string name)
{
    pFunc func(PyDict_GetItemString(main_dict, name.c_str()));
    return func;
}
