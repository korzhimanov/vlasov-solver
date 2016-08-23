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
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_sf_synchrotron.h>

/**
 * /class PFC
 * The class implements vlasov solver based on PFC scheme
*/
class PFC
{
private:
    double *f1,
           *f2;
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
    const double onesixth; // just a constant equal to 1/6
    double res, tmp, fold;
    int new_cell;
    double flux_part;
    double quart_q2n0dtdp_m, half_qn0dz, q2dt_dzdp, qdt_mdp, halfq_m, q_m, twicepisqrt3q2n0dpdzr_l, q2_m2, dp2, *p, *p2;
    pFunc *profile; // concentration profile function

public:
//------constructors---------------------------------------------------
    PFC();
    PFC(int particle_type, pyinput *in, Mesh *grid, double *n0, double *p0);
//------destructor-----------------------------------------------------
    ~PFC();

//------initializing---------------------------------------------------
    // initializes all parameters, allocates memory and fills it with zeros
    void Init(int particle_type, pyinput *in, Mesh *grid, double *n0, double *p0);
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
    void Destination(int old_cell, int &new_cell, int max_cell, double &flux_part);
    double SlopeCorrectorLeft(double centralvalue, double leftvalue);
    double SlopeCorrectorRight(double centralvalue, double rightvalue);

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
    void SaveEmitRadiationTxt(std::string filename, int max_harmonic, int dn, double dt, FDTD *fdtd, double *ez, double *ax, double *ay, double *a2);
    void SaveEmitRadiationBin(std::string filename, int max_harmonic, int dn, double dt, FDTD *fdtd, double *ez, double *ax, double *ay, double *a2);
    void SaveEmitRadiationGZip(std::string filename, int max_harmonic, int dn, double dt, FDTD *fdtd, double *ez, double *ax, double *ay, double *a2);
    // calculates kinetic energy of particles
    double KineticEnergy(double* ax, double *ay);
    double KineticEnergy(double* ax, double *ay, int left_bound, int right_bound);

private:
//------miscellaneous--------------------------------------------------
    void SetNotNegative(pyinput *in, std::string name, int *var);
    void SetPositive(pyinput *in, std::string name, int *var);
    void SetPositive(pyinput *in, std::string name, double *var);
    void SetFunction(pyinput *in, std::string name, pFunc &F);
    void CalcConcDstr(float* n);
    double CalcEmitRadiation(double freq, FDTD *fdtd, double *ez, double *ax, double *ay, double *a2);
};

inline void PFC::Destination(int old_cell, int &new_cell, int max_cell, double &flux_part)
{
    new_cell = old_cell;
    flux_part = 0.;

    if (vdt_dz > 0.)
    {
        new_cell = old_cell - 1 + (int)ceil(vdt_dz);
        flux_part = vdt_dz - floor(vdt_dz);
        if (new_cell >= max_cell - 1)
        {
            new_cell = max_cell - 1;
            flux_part = 0.;
        }
    }
    else
    {
        new_cell = old_cell + (int)ceil(vdt_dz);
        flux_part = vdt_dz - ceil(vdt_dz);
        if (new_cell <= 0)
        {
            new_cell = 0;
            flux_part = 0.;
        }
    }
}

// slope corrector for right boundary
inline double PFC::SlopeCorrectorRight(double centralvalue, double rightvalue)
{
    double res, tmp;
    if (fabs(rightvalue - centralvalue) < 1e-14)
        res = 0.;
    else
    {
        tmp = 1./(rightvalue - centralvalue);
        if (tmp > 0.)
        {
            if (2.*(centralvalue)*tmp < 1.)
                res = 2.*(centralvalue)*tmp;
            else
                res = 1.;
        }
        else
        {
            if (-2.*(maxf-centralvalue)*tmp < 1.)
                res = -2.*(maxf-centralvalue)*tmp;
            else
                res = 1.;
        }
    }
    return res;
}

// slope corrector for left boundary
inline double PFC::SlopeCorrectorLeft(double centralvalue, double leftvalue)
{
    double res, tmp;
    if (fabs(centralvalue - leftvalue) < 1e-14)
        res = 0.;
    else
    {
        tmp = 1./(centralvalue - leftvalue);
        if (tmp < 0.)
        {
            if (-2.*(centralvalue)*tmp < 1.)
                res = -2.*(centralvalue)*tmp;
            else
                res = 1.;
        }
        else
        {
            if (2.*(maxf-centralvalue)*tmp < 1.)
                res = 2.*(maxf-centralvalue)*tmp;
            else
                res = 1.;
        }
    }
    return res;
}

#endif // PFC_H
