/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file solver.h
 * \brief The header file which defines Solver class and implements its methods
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 * \bug The methods of the Solver class are implemented in the header due to some problems with compilation.
 */

#ifndef SOLVER_H
#define SOLVER_H

#ifdef _OPENMP
#include <omp.h>
#endif
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <math.h>
#include <cstdio>
#include "mymath.h"
#include "file_saving.h"
#include "initparams.h"
#include "fdtd.h"
#include "pfc.h"
#include "particle.h"
#include <string.h>

/**
 * \class Solver
 * The main class where actually the main sub-steps of each iteration are defined
 */
class Solver
{
    public:
        // initial parameters
        InitParams *params;

        // variables containing data
        Particle *prt; // particles
        double *fixed_ions_conc; // concentration of fixed ions
        FDTD *fdtd; // FDTD - class
        PFC *pfc; // PFC - class
        double *ax, *ay, *a2; // vector potential
        double *ez; // longitudinal field

        // auxiliary variables
        double em_flux; // flux density of electromagnetic field
        double electrostatic_energy, laser_energy, *plasma_energy, plasma_energy_0;

    public:
//------constructors---------------------------------------------------
        Solver();
        Solver(InitParams*);
//------destructor-----------------------------------------------------
        ~Solver();

//------initialisation-------------------------------------------------
        int InitVars(std::string, std::string);
        void AllocMemory();
        void InitFields();
        void InitPlasma(); // initializes plasma distribution
        void InitTestParticles();
        void CreateDirs();
        void SaveInput(std::string file);

//------calculations---------------------------------------------------
        void MoveParticles();
        inline void CalcSources(double t)
        {
            fdtd->exl =   params->pulse_x->call((double)(t - .5*params->mesh->dz));
            fdtd->eyl =   params->pulse_y->call((double)(t - .5*params->mesh->dz));
            fdtd->hxl = - params->pulse_y->call((double)(t + .5*params->mesh->dt));
            fdtd->hyl =   params->pulse_x->call((double)(t + .5*params->mesh->dt));
        }

        /**
         * \todo Change 0.5*pfc->N_0*params->mesh->dz to a single constant
         */
        inline void CalcLongFields()
        {
            // calculates longitudinal field
            memset(ez, 0, sizeof(double)*(params->mesh->MAX_Z+1));
            for (int sp = 0; sp < params->NUM_SP; sp++)
                pfc[sp].CalcLongitudinalField(ez);

            ez[0] += 0.5*pfc->N_0*params->mesh->dz*fixed_ions_conc[0];
            for (int i = 1; i < params->mesh->MAX_Z; i++)
                ez[i] += ez[i-1] + 0.5*pfc->N_0*params->mesh->dz*(fixed_ions_conc[i]+fixed_ions_conc[i-1]);
            ez[params->mesh->MAX_Z] += 0.5*pfc->N_0*params->mesh->dz*fixed_ions_conc[params->mesh->MAX_Z-1];
        }

        inline void CalcTransFields()
        {
            // evaluates vector potential
            for (int i = 0; i <= params->mesh->MAX_Z; i++)
            {
                ax[i] -= fdtd->ex[i]*params->mesh->dt;
                ay[i] -= fdtd->ey[i]*params->mesh->dt;
                a2[i] = ax[i]*ax[i]+ay[i]*ay[i];
            }

            // evaluating electric and magnetic fields
            fdtd->Maxwell();

            em_flux += fdtd->FluxIn();
        }

        /**
         * \todo Change 0.5*sin(params->THETA)*pfc->N_0*params->mesh->dz to a single constant
         */
        inline void FieldGeneration()
        {
            // field generation by plasma currents
            for (int sp = 0; sp < params->NUM_SP; sp++)
                pfc[sp].CalcCurrent(fdtd, ax, ay, a2);

            fdtd->ey[0] -= 0.5*sin(params->THETA)*pfc->N_0*params->mesh->dz*fixed_ions_conc[0];
            for (int i = 1; i < params->mesh->MAX_Z; i++)
                fdtd->ey[i] -= 0.5*sin(params->THETA)*pfc->N_0*params->mesh->dz*(fixed_ions_conc[i]+fixed_ions_conc[i-1]);
            fdtd->ey[params->mesh->MAX_Z] -= 0.5*sin(params->THETA)*pfc->N_0*params->mesh->dz*fixed_ions_conc[params->mesh->MAX_Z-1];
        }

        inline void CalcDstrFunc()
        {
            for (int sp = 0; sp < params->NUM_SP; sp++)
                pfc[sp].MakeStep(ez, ax, ay);
        }

//------saving data----------------------------------------------------
        void InitOutput(std::string);
        void SaveFields(int k); // saving fields data in files
        void SaveConcs(int k); // saving concentarations data in files
        void SaveDstrFunc(int k); // saving distribution functions data in files
        void SavePrtDat(int k); // save particles data
        void SaveOutput(int k, std::string file);
        void SaveResults();

    private:
};

Solver::Solver(InitParams *p)
{
    params = p;
    pfc = new PFC[params->NUM_SP];
    for (int sp = 0; sp < params->NUM_SP; sp++)
        pfc[sp].Init(sp, params->in, params->mesh);
    AllocMemory();
}

Solver::~Solver()
{
    fclose(params->output);
    delete[] pfc;

    delete[] fixed_ions_conc;

    delete fdtd;
    delete[] ax;
    delete[] ay;
    delete[] a2;
    delete[] ez;

    if (params->NUM_PRT > 0) delete[] prt;

    delete[] plasma_energy;
}

int Solver::InitVars(std::string file_name, std::string directory_name)
{
    // moved to InitParams::Init()
    return 0;
}

void Solver::AllocMemory()
{
    // memory allocation for fixed ions
    fixed_ions_conc = new double[params->mesh->MAX_Z];
    for (int sp = 0; sp < params->NUM_SP; sp++)
        pfc[sp].AllocMemory();

    // memory allocation for electromagnetic fields
    fdtd = new FDTD;
    ax = new double[params->mesh->MAX_Z + 1];
    ay = new double[params->mesh->MAX_Z + 1];
    a2 = new double[params->mesh->MAX_Z + 1];
    ez = new double[params->mesh->MAX_Z + 1];

    // memory allocation for test particles
    if (params->NUM_PRT > 0) prt = new Particle[params->NUM_PRT];

    plasma_energy = new double[params->NUM_SP];
}

void Solver::InitOutput(std::string fn)
{
    params->output = fopen((params->output_directory_name + fn).c_str(), "w+");
    fprintf(params->output, " Step |   Neutrality   |");
    for (int sp = 0; sp < params->NUM_SP; sp++)
        fprintf(params->output, "       W_%d       |", sp);
    fprintf(params->output, "     W_el-st    |     W_laser    |     EM Flux    |     W_total    |\n");
    fflush(params->output);
}

void Solver::InitFields()
{
    fdtd->Init(params->mesh, params->source);
    fdtd->InitPML(int(1./params->mesh->dz), 10.);
    mymath::zeros(ax, params->mesh->MAX_Z+1);
    mymath::zeros(ay, params->mesh->MAX_Z+1);
    mymath::zeros(a2, params->mesh->MAX_Z+1);
    mymath::zeros(ez, params->mesh->MAX_Z+1);

    em_flux = 0.;
}

void Solver::InitPlasma()
{
    for (int i = 0; i < params->mesh->MAX_Z; i++)
        fixed_ions_conc[i] = params->fixed_ions_profile->call(i);
    for (int sp = 0; sp < params->NUM_SP; sp++)
        pfc[sp].SetDistribution();
}

void Solver::InitTestParticles()
{
    if (params->NUM_PRT > 0) for (int i = 0; i < params->NUM_PRT; i++)
    {
        prt[i].q = int(params->CHARGE_PRT);
        prt[i].m = params->MASS_PRT;
        for (int j = 0; j < 3; j++)
        {
            prt[i].r[j] = 0.;
            prt[i].p[j] = 0.;
        }
        prt[i].p[1] = tan(params->THETA);
        prt[i].r[2] = (params->start_point + params->interval*i)*params->mesh->dz;
    }
}

void Solver::CreateDirs()
{
    // creating main directory for output
    {
        filesaving::create_dir("", (params->output_directory_name).c_str());
    }

    // creating directory for test particles
    if (params->NUM_PRT > 0)
    {
        filesaving::create_dir(params->output_directory_name, "particles/");
    }

    // creating directories for storing distribution functions in files
    if (params->save_dstr == true)
        for (int sp = 0; sp < params->NUM_SP; sp++)
        {
            filesaving::create_dir(params->output_directory_name, "dstr_func%d/", sp);
        }
}

void Solver::SaveInput(std::string file)
{
    // storing information about input parameters
    FILE *input;
    input = filesaving::open_file("w+", params->output_directory_name, file.c_str());
    fprintf(input, "\nVlasov-Maxwell code. Developed by Artem Korzhimanov and Arkady Gonoskov. mailto: kav@ufp.appl.sci-nnov.ru\n");

    fprintf(input, "\nGENERAL PARAMETERS\n\n");
    fprintf(input, "\tSpace cells per wavelength = %d\n", params->ppw);
    fprintf(input, "\tSpace step = %f Wavelengths\n", .5*M_1_PI*params->mesh->dz);
    fprintf(input, "\tTotal space cells = %d\n", params->mesh->MAX_Z);
    fprintf(input, "\tTotal space size = %f Wavelengths\n", (double)params->mesh->MAX_Z/params->ppw);
    fprintf(input, "\tTime step = %f Waveperiods\n", .5*M_1_PI*params->mesh->dt);
    fprintf(input, "\tTotal time steps = %d\n", params->mesh->MAX_T);
    fprintf(input, "\tTotal running time = %f Waveperiods\n", .5*M_1_PI*params->mesh->MAX_T*params->mesh->dt);

    fprintf(input, "\nPLASMA PARAMETERS\n\n");
    for (int sp = 0; sp < params->NUM_SP; sp++)
    {
        fprintf(input, "\n\tSpecies %d\n\n", sp);
        pfc[sp].SaveInput(input);
    }

    fprintf(input, "\nTEST PARTICLES PARAMETERS\n\n");
    fprintf(input, "\tNumber of particles = %d\n", params->NUM_PRT);
    if (params->NUM_PRT > 0)
    {
        fprintf(input, "\tMass of particles = %f\n", params->MASS_PRT);
        fprintf(input, "\tCharge of particles = %f\n", params->CHARGE_PRT);
        fprintf(input, "\tParticles are equidistantly deposited from %f to %f Wavelenths\n", (float)params->start_point/params->ppw, (float)(params->start_point + params->interval*params->NUM_PRT)/params->ppw);
    }

    fprintf(input, "\nOUTPUT PARAMETERS\n\n");
    if (params->save_fields == true)
        fprintf(input, "\tFields are being saved every %d steps or every %f Waveperiods in files in %s-format\n", params->save_fields_dt, (double)params->save_fields_dt/params->ppw, params->save_fields_format.c_str());
    else fprintf(input, "\tFields are not being saved\n");
    if (params->save_concs == true)
        fprintf(input, "\tConcentration distributions are being saved every %d steps or every %f Waveperiods in files in %s-format\n", params->save_concs_dt, (double)params->save_concs_dt/params->ppw, params->save_concs_format.c_str());
    else fprintf(input, "\tConcentration distributions are not being saved\n");
    if (params->save_dstr == true)
        fprintf(input, "\tDistribution functions are being saved every %d steps or every %f Waveperiods in files in %s-format\n", params->save_dstr_dt, (double)params->save_dstr_dt/params->ppw, params->save_dstr_format.c_str());
    else fprintf(input, "\tDistribution functions are not being saved\n");

    filesaving::close_file(input);
}

void Solver::MoveParticles()
{
    int z;
    double dz;

    int i;

    double *E = new double[3];
    double *B = new double[3];
    B[2] = 0;
    for (i = 0; i < params->NUM_PRT; i++)
    {
        z = int(floor(prt[i].r[2]/params->mesh->dz));
        if (z < 0 || z >= params->mesh->MAX_Z) continue;
        dz = prt[i].r[2]/params->mesh->dz - floor(prt[i].r[2]/params->mesh->dz);

        E[0] = (1. - dz) * fdtd->ex[z] + dz * fdtd->ex[z+1];
        E[1] = (1. - dz) * fdtd->ey[z] + dz * fdtd->ey[z+1];

        z = int(floor(prt[i].r[2]/params->mesh->dz - 0.5));
        if (z < 0 || z >= params->mesh->MAX_Z-1) continue;
        dz = prt[i].r[2]/params->mesh->dz - floor(prt[i].r[2]/params->mesh->dz);

        E[2] = (1. - dz) * ez[z] + dz * ez[z+1];
        B[0] = (1. - dz) * fdtd->hx[z] + dz * fdtd->hx[z+1];
        B[1] = (1. - dz) * fdtd->hy[z] + dz * fdtd->hy[z+1];

        prt[i].MakeStep(params->mesh->dt, E, B);
    }

    delete[] E;
    delete[] B;
}

void Solver::SaveFields(int k)
{
    if (params->save_fields == true && k%params->save_fields_dt == 0)
    {
        std::string format;
        if (params->save_fields_format == "")
            format = params->save_format;
        else
            format = params->save_fields_format;

        if (format == "txt")
        {
            filesaving::save_file_1D(fdtd->ex, params->mesh->MAX_Z, params->output_directory_name, "ex.txt");
            filesaving::save_file_1D(fdtd->ey, params->mesh->MAX_Z, params->output_directory_name, "ey.txt");
            filesaving::save_file_1D(fdtd->hx, params->mesh->MAX_Z, params->output_directory_name, "hx.txt");
            filesaving::save_file_1D(fdtd->hy, params->mesh->MAX_Z, params->output_directory_name, "hy.txt");
            filesaving::save_file_1D(ez, params->mesh->MAX_Z, params->output_directory_name, "ez.txt");

            filesaving::save_file_1D(a2, params->mesh->MAX_Z, params->output_directory_name, "vecpot2.txt");
        }

        if (format == "bin")
        {
            filesaving::save_file_1D_bin(fdtd->ex, params->mesh->MAX_Z, params->output_directory_name, "ex.bin");
            filesaving::save_file_1D_bin(fdtd->ey, params->mesh->MAX_Z, params->output_directory_name, "ey.bin");
            filesaving::save_file_1D_bin(fdtd->hx, params->mesh->MAX_Z, params->output_directory_name, "hx.bin");
            filesaving::save_file_1D_bin(fdtd->hy, params->mesh->MAX_Z, params->output_directory_name, "hy.bin");
            filesaving::save_file_1D_bin(ez, params->mesh->MAX_Z, params->output_directory_name, "ez.bin");

            filesaving::save_file_1D_bin(a2, params->mesh->MAX_Z, params->output_directory_name, "vecpot2.bin");
        }

        if (format == "gzip")
        {
            filesaving::save_file_1D_gzip(fdtd->ex, params->mesh->MAX_Z, params->output_directory_name, "ex.gz");
            filesaving::save_file_1D_gzip(fdtd->ey, params->mesh->MAX_Z, params->output_directory_name, "ey.gz");
            filesaving::save_file_1D_gzip(fdtd->hx, params->mesh->MAX_Z, params->output_directory_name, "hx.gz");
            filesaving::save_file_1D_gzip(fdtd->hy, params->mesh->MAX_Z, params->output_directory_name, "hy.gz");
            filesaving::save_file_1D_gzip(ez, params->mesh->MAX_Z, params->output_directory_name, "ez.gz");

            filesaving::save_file_1D_gzip(a2, params->mesh->MAX_Z, params->output_directory_name, "vecpot2.gz");
        }

    }
}

void Solver::SaveConcs(int k)
{
    if (params->save_concs == true && k%params->save_concs_dt == 0)
    {
        std::string format;
        if (params->save_concs_format == "")
            format = params->save_format;
        else
            format = params->save_concs_format;

        for (int sp = 0; sp < params->NUM_SP; sp++)
        {
            std::stringstream ss;
            ss << params->output_directory_name << "conc" << sp;
            if (format == "txt")
                pfc[sp].SaveConcentrationTxt(ss.str());

            if (format == "bin")
                pfc[sp].SaveConcentrationBin(ss.str());

            if (format == "gzip")
                pfc[sp].SaveConcentrationGZip(ss.str());
        }
    }
}

void Solver::SaveDstrFunc(int k)
{
    if (params->save_dstr == true && k%params->save_dstr_dt == 0)
    {
        std::string format;
        if (params->save_dstr_format == "")
            format = params->save_format;
        else
            format = params->save_dstr_format;

        for (int sp = 0; sp < params->NUM_SP; sp++)
        {
            std::stringstream ss;
            ss << params->output_directory_name << "dstr_func" << sp << "/data" << std::setfill('0') << std::setw(8) << k;
            if (format == "txt")
                pfc[sp].SaveDstrFunctionTxt(ss.str());

            if (format == "bin")
                pfc[sp].SaveDstrFunctionBin(ss.str());

            if (format == "gzip")
                pfc[sp].SaveDstrFunctionGZip(ss.str());
        }
    }
}

void Solver::SavePrtDat(int k)
{
    if (params->NUM_PRT > 0)
    {
        std::ofstream *out_file;
        for (int i = 0; i < params->NUM_PRT; i++)
        {
            char file_path[256];
            strcpy(file_path, params->output_directory_name.c_str());
            char file_name[256];
            sprintf(file_name, "particles/particle%06d.txt", i);
            strcat(file_path, file_name);
            out_file = new std::ofstream(file_path, std::ios_base::app);
            *out_file << prt[i].r[0] << "\t" << prt[i].r[1] << "\t" << prt[i].r[2] << "\t" << prt[i].p[0] << "\t" << prt[i].p[1] << "\t" << prt[i].p[2] << std::endl;
            out_file->close();
            delete out_file;
        }
    }
}

void Solver::SaveOutput(int k, std::string file)
{
    if (k%params->ppw == 0)
    {
        for (int sp = 0; sp < params->NUM_SP; sp++)
            plasma_energy[sp] = pfc[sp].KineticEnergy(ax, ay);

        laser_energy = fdtd->Energy();

        electrostatic_energy = 0.;
        for (int i = 0; i <= params->mesh->MAX_Z; i++)
            electrostatic_energy += ez[i]*ez[i];
        electrostatic_energy *= 0.5*params->mesh->dz;

        FILE *output;
        output = filesaving::open_file("a+", params->output_directory_name, file.c_str());
        fprintf(output, "%06d|%16.8e|", k, ez[params->mesh->MAX_Z-1]);
        for (int sp = 0; sp < params->NUM_SP; sp++)
            fprintf(output, " %16.8e|", plasma_energy[sp]);
        fprintf(output, "%16.8e|%16.8e|%16.8e|%16.8e|\n", electrostatic_energy, laser_energy, em_flux, mymath::sum(plasma_energy, params->NUM_SP)+electrostatic_energy+laser_energy-em_flux);
        filesaving::close_file(output);

        if (k==0) plasma_energy_0 = mymath::sum(plasma_energy, params->NUM_SP);
    }
}

void Solver::SaveResults()
{
    FILE *en_out;
    en_out = filesaving::open_file("a+", "", "energy.txt");
    double plasma_energy_1 = mymath::sum(plasma_energy, params->NUM_SP)+electrostatic_energy+laser_energy;
    fprintf(en_out, "%.8g\t%.8g\t%.8g\n", plasma_energy_0, plasma_energy_1, 2.*(plasma_energy_1-plasma_energy_0));
    filesaving::close_file(en_out);

    en_out = filesaving::open_file("a+", "", "energy2.txt");
    fprintf(en_out, "%.8g\t", (plasma_energy_1-plasma_energy_0)/(0.5));
    filesaving::close_file(en_out);
}

#endif // SOLVER_H
