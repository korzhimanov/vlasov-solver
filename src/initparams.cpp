/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file initparams.cpp
 * \brief The file implements methods of InitParams class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 *
 * The file contains argument parser, timer, main cycle and basic logging.
 */


#include "initparams.h"
#include <string>
#include <iostream>
#include "file_saving.h"

InitParams::InitParams(pyinput* i, std::string dn, Mesh* m, int* err) : mesh(m), in(i)
{
    *err = Init(dn);
}

InitParams::~InitParams()
{
}

int InitParams::Init(std::string dn)
{
    output_directory_name = dn + "/";

    int tpi;
    double tpf;

    THETA = in->GetDouble("THETA");
    if (THETA < 0. && THETA >= 0.5*M_PI)
        std::cout << "Invalid incident angle! It should be between 0 and pi/2." << std::endl;

    tpi = 0; if ( !in->SetNotNegative("NUM_PRT", &tpi) ) return  1600; NUM_PRT = tpi;
    if (NUM_PRT > 0)
    {
        tpi =  0; if ( !in->SetNotNegative("start_point", &tpi) ) return  1610; start_point = tpi;
        tpf    = 1.; if ( !in->SetPositive   ("interval", &tpf) ) return  1620; interval = tpf;
        tpf    = 1.; if ( !in->SetPositive   ("MASS_PRT", &tpf) ) return  1630; MASS_PRT = tpf;
        CHARGE_PRT = in->GetDouble("CHARGE_PRT");
    }

    save_fields = in->GetInt("save_fields");
    save_concs  = in->GetInt("save_concs");
    save_dstr   = in->GetInt("save_dstr");

    if ( (save_fields || save_concs || save_dstr) == true)
    {
        tpi = mesh->ppw; if ( !in->SetPositive("save_dt", &tpi) ) return 1410; save_dt = tpi;
        save_format = "txt"; if (   !SetFormat("save_format", &(save_format))) return 1420;

        if (save_fields == true)
        {
            tpi = 0; if ( !in->SetNotNegative("save_fields_dt", &tpi) ) return 1430; save_fields_dt = tpi;
            if (save_fields_dt == 0)
                save_fields_dt = save_dt;
            save_fields_format = ""; if ( !SetFormat("save_fields_format", &(save_fields_format)) ) return 1431;
            if (save_fields_format == "")
                save_fields_format = save_format;
        }
        if (save_concs == true)
        {
            tpi = 0; if ( !in->SetNotNegative(    "save_concs_dt", &tpi)    ) return 1440; save_concs_dt = tpi;
            if (save_concs_dt == 0)
                save_concs_dt = save_dt;
            save_concs_format = ""; if ( !SetFormat("save_concs_format", &(save_concs_format)) ) return 1441;
            if (save_concs_format == "")
                save_concs_format = save_format;
        }
        if (save_dstr == true)
        {
            tpi = 0; if ( !in->SetNotNegative("save_dstr_dt", &tpi) ) return 1450; save_dstr_dt = tpi;
            if (save_dstr_dt == 0)
                save_dstr_dt = save_dt;
            save_dstr_format = ""; if ( !SetFormat("save_dstr_format", &(save_dstr_format)) ) return 1451;
            if (save_dstr_format == "")
                save_dstr_format = save_format;
        }
    }

    return 0;
}

bool InitParams::SetFormat(std::string name, std::string *var)
{
    *var = in->GetString(name);
    if (save_format == "txt" || save_format == "bin" || save_format ==  "gzip") return true;
    else {std::cout << name + " must be either 'txt' or 'bin' or 'gzip'" << std::endl; return false;}
}