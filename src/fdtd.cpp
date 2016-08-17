/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file fdtd.cpp
 * \brief The source file which implements the methods of the FDTD class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#ifdef _OPENMP
#include <omp.h>
#endif
#include <iostream>
#include "fdtd.h"

FDTD::FDTD()
{
    exl = eyl = exr = eyr = hxl = hyl = hxr = hyr = 0.;
}

FDTD::FDTD(Mesh *grid)
{
    exl = eyl = exr = eyr = hxl = hyl = hxr = hyr = 0.;
    Init(grid);
    InitPML(512, 100.);
}

FDTD::~FDTD()
{
    delete[] ex;
    delete[] ey;
    delete[] hx;
    delete[] hy;
    delete[] r;
    delete[] r1;
}

void FDTD::Init(Mesh *grid, int source)
{
    mesh = grid;
    dtdz = mesh->dt/mesh->dz;
    ex = new double[mesh->MAX_Z + 1];
    ey = new double[mesh->MAX_Z + 1];
    hx = new double[mesh->MAX_Z];
    hy = new double[mesh->MAX_Z];
    for (int i = 0; i < mesh->MAX_Z; i++)
        ex[i] = ey[i] = hx[i] = hy[i] = 0.;
    ex[mesh->MAX_Z] = ey[mesh->MAX_Z] = 0.;
    SOURCE = source;
}

void FDTD::InitPML(int pml, double max_sigma)
{
    PML = pml;
    MAX_SIGMA = max_sigma;
    r  = new double[2*PML];
    r1 = new double[2*PML];
    CalcPMLCoeff();
}

void FDTD::Maxwell()
{
    // common Yee's algorithm
    int j, p;

    // advance magnetic fields
    for (j = PML; j <= mesh->MAX_Z - PML; j++)
    {
        hx[j] += dtdz * (ey[j+1] - ey[j]);
        hy[j] -= dtdz * (ex[j+1] - ex[j]);
    }

    // PML with exponential time-stepping
    for (j = 0, p = 2*PML - 1; j <= PML - 1; j++, p-=2)
    {
        hx[j] = r[p] * hx[j] - r1[p] * (ey[j+1] - ey[j]);
        hy[j] = r[p] * hy[j] + r1[p] * (ex[j+1] - ex[j]);
    }

    for (j = mesh->MAX_Z - PML + 1, p = 1; j < mesh->MAX_Z; j++, p+=2)
    {
        hx[j] = r[p] * hx[j] - r1[p] * (ey[j+1] - ey[j]);
        hy[j] = r[p] * hy[j] + r1[p] * (ex[j+1] - ex[j]);
    }

    // incident wave correction
    hx[PML+SOURCE-1]      -= dtdz * eyl;
    hy[PML+SOURCE-1]      += dtdz * exl;
    hx[mesh->MAX_Z-PML-1] += dtdz * eyr;
    hy[mesh->MAX_Z-PML-1] -= dtdz * exr;

    // advance electric fields
    for (j = PML + 1; j <= mesh->MAX_Z - PML; j++)
    {
        ex[j] -= dtdz * (hy[j] - hy[j-1]);
        ey[j] += dtdz * (hx[j] - hx[j-1]);
    }

    // PML with exponential time-stepping
    for (j = 1, p = 2*PML-2; j <= PML; j++, p-=2)
    {
        ex[j] = r[p] * ex[j] + r1[p] * (hy[j] - hy[j-1]);
        ey[j] = r[p] * ey[j] - r1[p] * (hx[j] - hx[j-1]);
    }

    for (j = mesh->MAX_Z - PML + 1, p = 0; j < mesh->MAX_Z; j++, p+=2)
    {
        ex[j] = r[p] * ex[j] + r1[p] * (hy[j] - hy[j-1]);
        ey[j] = r[p] * ey[j] - r1[p] * (hx[j] - hx[j-1]);
    }

    // incident wave correction
    ex[PML+SOURCE]    += dtdz * hyl;
    ey[PML+SOURCE]    -= dtdz * hxl;
    ex[mesh->MAX_Z-PML] += dtdz * hyr;
    ey[mesh->MAX_Z-PML] -= dtdz * hxr;
}

void FDTD::CalcPMLCoeff()
{
    double sigma;
    for (int j = 0; j < 2*PML; j++)
    {
        sigma = MAX_SIGMA * pow(double(j+1) / double(2*PML), 2.);
        r[j] = exp( - sigma * mesh->dt);
        r1[j] = (r[j] - 1.) / (sigma * mesh->dz);
    }
}
