/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file pfc.h
 * \brief The header file which defines PFC class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#ifndef PFC_H
#define PFC_H

#ifdef _OPENMP
#include <omp.h>
#endif
#include <cstdio>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <zlib.h>
#include "pyinput.h"
#include "mesh.h"
#include "fdtd.h"
#include "mymath.h"
#include "file_saving.h"

/**
 * /class PFC
 * The class implements vlasov solver based on PFC scheme
*/
class PFC
{
public:
//------constructors---------------------------------------------------
    PFC();
    PFC(int particle_type, pyinput *in, Mesh *grid, double *n0, double *p0, int *err);
//------destructor-----------------------------------------------------
    ~PFC();

//------initializing---------------------------------------------------
    // initializes all parameters, allocates memory and fills it with zeros
    int Init(int particle_type, pyinput *in, Mesh *grid, double *n0, double *p0);
    void AllocMemory();
    void SetDistribution();

//------solver---------------------------------------------------------
    // makes step
    void MakeStep(double *ez, double *ax, double *ay);
    // calculates longitudinal field
    void CalcLongitudinalField(double *ez);
    void CalcCurrent(FDTD *fdtd, double *ax, double *ay, double *a2);

private:
    // determines destination cell and parts of flux between it and nighbouring one
    double CorrectedLeftSlope (double *centralvalue, double *leftvalue);
    double CorrectedRightSlope(double *centralvalue, double *rightvalue);

//------informative functions------------------------------------------
public:
    // saves input information
    void SaveInput(std::ofstream&);
    // saves output information
    void SaveConcentrationTxt(std::string filename);
    void SaveConcentrationBin(std::string filename);
    void SaveConcentrationGZip(std::string filename);
    void SaveDstrFunctionTxt(std::string filename);
    void SaveDstrFunctionBin(std::string filename);
    void SaveDstrFunctionGZip(std::string filename);
    // calculates kinetic energy of particles
    double KineticEnergy(double* ax, double *ay);
    double KineticEnergy(double* ax, double *ay, int left_bound, int right_bound);

private:
//------miscellaneous--------------------------------------------------
    void CalcConcDstr();

private:
    double *f1,
           *f2,
           *conc;
    Mesh *mesh;
    double N0, // critical concentration in a boosted frame (useful for simulation of oblique irradiation)
           P0; // initial mean momentum along y-axis (useful for simulation of oblique irradiation in boosted frame)
    int MAX_P, // number of steps in momentum space
        type; // type of the particle
    double
        MASS, // particle mass
        CHARGE, // particle charge
        T_init, // initial temperature
        MEAN_P, // mean momentum along z-axis
        dp; // mesh step in momentum space
    double maxf; // maximal value of the distribution function
    double vdt_dz; // projections of velocity in phase space
    double fold, flux_part;
    int new_cell;
    double n0dp, q_m, q2_m2, halfq_m, halfqdz, quart_q2n0dtdp_m, qdt_mdp, q2dt_dzdp, *p, *p2;
    pFunc *profile; // concentration profile function
};
#endif // PFC_H
