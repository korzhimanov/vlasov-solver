/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file initparams.h
 * \brief The file defines InitParams class and Params structure
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 *
 * The file contains argument parser, timer, main cycle and basic logging.
 */


#ifndef INITPARAMS_H
#define INITPARAMS_H

#include <string>
#include "pfunc.h"
#include "pyinput.h"
#include "mesh.h"

class InitParams
{
    public:
        // general
        Mesh *mesh;
        double THETA;
        // particles
        int NUM_PRT,
            start_point;
        double MASS_PRT,
               CHARGE_PRT,
               interval;
        // output
        std::string output_directory_name;
        FILE *output;
        bool save_fields,
             save_concs,
             save_dstr;
        int save_dt,
            save_fields_dt,
            save_concs_dt,
            save_dstr_dt;
        std::string save_format,
                    save_fields_format,
                    save_concs_format,
                    save_dstr_format;
        pyinput *in; // move to private later

    public:
        InitParams(pyinput*, std::string, Mesh*, int*);
        virtual ~InitParams();

    private:
        int Init(std::string);
        bool SetFormat(std::string name, std::string *var);
};

#endif // INITPARAMS_H
