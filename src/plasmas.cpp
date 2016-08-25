/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file plasmas.cpp
 * \brief The source file which implements the methods of the Plasmas class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */


#include "plasmas.h"
#include "mymath.h"
#include "pfunc.h"
#include "pyinput.h"
#include "pfc.h"

Plasmas::Plasmas(pyinput* in, Mesh* m, int* err) : species_number(0), critical_concentration(1.), mean_initial_momentum(0.), mesh(m)
{
    *err = Init(in);
    AllocMemory();
}

Plasmas::~Plasmas()
{
    if (fixed_ions_profile) delete fixed_ions_profile;
    delete[] fixed_ions_conc;
    delete[] pfc;
    delete[] ax;
    delete[] ay;
    delete[] a2;
}

int Plasmas::Init(pyinput* in)
{
    fixed_ions_profile = new pFunc(in->GetFunc("FIXED_IONS_PROFILE"));

    if ( !in->SetNotNegative("NUM_SP", &species_number) ) return 200;

    double THETA = 0.;
    THETA = in->GetDouble("THETA");
    if (THETA >= 0. && THETA < M_PI/2)
    {
        mean_initial_momentum = tan(THETA);
        critical_concentration = 1./(cos(THETA)*cos(THETA)*cos(THETA));
    }
    else return 210;

    pfc = new PFC[species_number];
    for (int sp = 0; sp < species_number; sp++)
    {
        int err = pfc[sp].Init(sp, in, mesh, &critical_concentration, &mean_initial_momentum);
        if (err != 0) return err;
    }

    return 0;
}

void Plasmas::AllocMemory()
{
    fixed_ions_conc = new double[mesh->MAX_Z];

    for (int sp = 0; sp < species_number; sp++)
        pfc[sp].AllocMemory();

    ax = new double[mesh->MAX_Z + 1];
    ay = new double[mesh->MAX_Z + 1];
    a2 = new double[mesh->MAX_Z + 1];
    ez = new double[mesh->MAX_Z + 1];
}

void Plasmas::InitDistribution()
{
    for (int i = 0; i < mesh->MAX_Z; i++)
        fixed_ions_conc[i] = fixed_ions_profile->call(i);
    delete fixed_ions_profile; // not needed anymore

    for (int sp = 0; sp < species_number; sp++)
        pfc[sp].SetDistribution();

    mymath::zeros(ax, mesh->MAX_Z+1);
    mymath::zeros(ay, mesh->MAX_Z+1);
    mymath::zeros(a2, mesh->MAX_Z+1);
    mymath::zeros(ez, mesh->MAX_Z+1);
}
