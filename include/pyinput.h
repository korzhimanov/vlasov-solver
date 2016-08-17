/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file pyinput.h
 * \brief The header file which defines pyinput class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#ifndef PYINPUT_H
#define PYINPUT_H

#include "pfunc.h"
#include <Python.h>
#include <iostream>

/**
 * /class pFunc
 * The class for reading python script with input parameters
*/
class pyinput
{
    public:
        pyinput();
        virtual ~pyinput();

        void ReadFile(std::string);

        int GetInt(std::string);
        double GetDouble(std::string);
        std::string GetString(std::string);

        pFunc GetFunc(std::string);

    private:
        FILE *f;
        PyObject *main_module, *main_dict;
};

#endif // PYINPUT_H
