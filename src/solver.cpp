/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file solver.cpp
 * \brief The main solver class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 * \bug The methods of the Solver class are implemented in the header due to some problems with compilation.
 */

#include "solver.h"
#include <cmath>
#include "pyinput.h"
#include "mesh.h"
#include "plasmas.h"
#include "fdtd.h"
#include "testparticles.h"

Solver::Solver(pyinput* in, Mesh* m, int* err) : mesh(m)
{
    plasmas = new Plasmas(in, m, err);
    fdtd = new FDTD(in, m, err);
    particles = new TestParticles(in, err);
    Init(in, err);
}

Solver::~Solver()
{
    delete plasmas;
    delete fdtd;
    delete particles;
}

int Solver::Init(pyinput* in, int* err)
{
    double THETA = 0.;
    THETA = in->GetDouble("THETA");
    if (THETA >= 0. && THETA < M_PI/2)
    {
        halfn0dz = 0.5*sin(THETA)*plasmas->critical_concentration*mesh->dz;
    }
    else return 10;

    return 0;
}

/**
 * \todo Remove z
 */
void Solver::MoveParticles()
{
    int z;
    double dz;

    int i;

    double *E = new double[3];
    double *B = new double[3];
    B[2] = 0;
    for (i = 0; i < particles->particles_number; i++)
    {
        z = int(floor(particles->prt[i].r[2]/mesh->dz));
        if (z < 0 || z >= mesh->MAX_Z) continue;
        dz = particles->prt[i].r[2]/mesh->dz - floor(particles->prt[i].r[2]/mesh->dz);

        E[0] = (1. - dz) * fdtd->ex[z] + dz * fdtd->ex[z+1];
        E[1] = (1. - dz) * fdtd->ey[z] + dz * fdtd->ey[z+1];

        z = int(floor(particles->prt[i].r[2]/mesh->dz - 0.5));
        if (z < 0 || z >= mesh->MAX_Z-1) continue;
        dz = particles->prt[i].r[2]/mesh->dz - floor(particles->prt[i].r[2]/mesh->dz);

        E[2] = (1. - dz) * plasmas->ez[z] + dz * plasmas->ez[z+1];
        B[0] = (1. - dz) * fdtd->hx[z] + dz * fdtd->hx[z+1];
        B[1] = (1. - dz) * fdtd->hy[z] + dz * fdtd->hy[z+1];

        particles->prt[i].MakeStep(mesh->dt, E, B);
    }

    delete[] E;
    delete[] B;
}
