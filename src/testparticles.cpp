/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file testparticles.cpp
 * \brief The source file which implements the methods of the TestParticles class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */


#include "testparticles.h"
#include <cmath>
#include "mesh.h"
#include "particle.h"
#include "pyinput.h"
#include "pfunc.h"


TestParticles::TestParticles(pyinput* in, int* err) : particles_number(0), start_point(0), interval(1.), mass(1.), charge(0.), mean_initial_momentum(0.)
{
    *err = Init(in);
}

TestParticles::~TestParticles()
{
    if (particles_number > 0) delete[] prt;
}

int TestParticles::Init(pyinput* in)
{
    if ( !in->SetNotNegative("NUM_PRT", &particles_number) ) return  400;
    if (particles_number > 0)
    {
        if ( !in->SetNotNegative("start_point", &start_point) ) return  410;
        if ( !in->SetPositive   ("interval", &interval) ) return  420;
        if ( !in->SetPositive   ("MASS_PRT", &mass) ) return  430;
        charge = in->GetDouble("CHARGE_PRT");
    }

    double THETA = 0.;
    THETA = in->GetDouble("THETA");
    if (THETA >= 0. && THETA < M_PI/2)
    {
        mean_initial_momentum = tan(THETA);
    }
    else return 440;

    return 0;
}

void TestParticles::AllocMemory()
{
    if (particles_number > 0) prt = new Particle[particles_number];
}

void TestParticles::InitParticles(Mesh* mesh)
{
    if (particles_number > 0) for (int i = 0; i < particles_number; i++)
    {
        prt[i].q = int(charge);
        prt[i].m = mass;
        for (int j = 0; j < 3; j++)
        {
            prt[i].r[j] = 0.;
            prt[i].p[j] = 0.;
        }
        prt[i].p[1] = mean_initial_momentum;
        prt[i].r[2] = (start_point + interval*i)*mesh->dz;
    }
}
