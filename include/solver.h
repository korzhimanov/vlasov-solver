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
#include "plasmas.h"
#include "particle.h"
#include "testparticles.h"
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
        Plasmas *plasmas; // Plasmas
        FDTD *fdtd; // FDTD - class
        TestParticles *particles;

        // auxiliary variables
        double em_flux; // flux density of electromagnetic field
        double electrostatic_energy, laser_energy, *plasma_energy, plasma_energy_0;

    public:
//------constructors---------------------------------------------------
        Solver(InitParams*, pyinput*, Mesh*, int*);
//------destructor-----------------------------------------------------
        ~Solver();

//------initialisation-------------------------------------------------
        int InitVars(std::string, std::string);
        void AllocMemory();
        void CreateDirs();
        void SaveInput(std::string file);

//------calculations---------------------------------------------------
        inline void CalcTransFields()
        {
            // evaluates vector potential
            for (int i = 0; i <= params->mesh->MAX_Z; i++)
            {
                plasmas->ax[i] -= fdtd->ex[i]*params->mesh->dt;
                plasmas->ay[i] -= fdtd->ey[i]*params->mesh->dt;
                plasmas->a2[i] = plasmas->ax[i]*plasmas->ax[i]+plasmas->ay[i]*plasmas->ay[i];
            }

            // evaluating electric and magnetic fields
            fdtd->Maxwell();

            em_flux += fdtd->FluxIn();
        }

        void MoveParticles();

        /**
         * \todo Change 0.5*sin(params->THETA)*pfc->N_0*params->mesh->dz to a single constant
         */
        inline void FieldGeneration()
        {
            // field generation by plasma currents
            for (int sp = 0; sp < plasmas->species_number; sp++)
                plasmas->pfc[sp].CalcCurrent(fdtd, plasmas->ax, plasmas->ay, plasmas->a2);

            // field generation by fixed ions (only in a boosted frame)
            fdtd->ey[0] -= 0.5*sin(params->THETA)*plasmas->critical_concentration*params->mesh->dz*plasmas->fixed_ions_conc[0];
            for (int i = 1; i < params->mesh->MAX_Z; i++)
                fdtd->ey[i] -= 0.5*sin(params->THETA)*plasmas->critical_concentration*params->mesh->dz*(plasmas->fixed_ions_conc[i]+plasmas->fixed_ions_conc[i-1]);
            fdtd->ey[params->mesh->MAX_Z] -= 0.5*sin(params->THETA)*plasmas->critical_concentration*params->mesh->dz*plasmas->fixed_ions_conc[params->mesh->MAX_Z-1];
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

Solver::Solver(InitParams* p, pyinput* in, Mesh* m, int* err) : params(p)
{
    plasmas = new Plasmas(in, m, err);
    fdtd = new FDTD(in, m, err);
    particles = new TestParticles(in, err);
    AllocMemory();
}

Solver::~Solver()
{
    fclose(params->output);

    delete plasmas;

    delete fdtd;

    delete particles;

    delete[] plasma_energy;
}

int Solver::InitVars(std::string file_name, std::string directory_name)
{
    // moved to InitParams::Init()
    return 0;
}

void Solver::AllocMemory()
{
    plasma_energy = new double[plasmas->species_number];
}

void Solver::InitOutput(std::string fn)
{
    em_flux = 0.;
    params->output = fopen((params->output_directory_name + fn).c_str(), "w+");
    fprintf(params->output, " Step |   Neutrality   |");
    for (int sp = 0; sp < plasmas->species_number; sp++)
        fprintf(params->output, "       W_%d       |", sp);
    fprintf(params->output, "     W_el-st    |     W_laser    |     EM Flux    |     W_total    |\n");
    fflush(params->output);
}

void Solver::CreateDirs()
{
    // creating main directory for output
    {
        filesaving::create_dir("", (params->output_directory_name).c_str());
    }

    // creating directory for test particles
    if (particles->particles_number > 0)
    {
        filesaving::create_dir(params->output_directory_name, "particles/");
    }

    // creating directories for storing distribution functions in files
    if (params->save_dstr == true)
        for (int sp = 0; sp < plasmas->species_number; sp++)
        {
            filesaving::create_dir(params->output_directory_name, "dstr_func%d/", sp);
        }
}

void Solver::SaveInput(std::string file)
{
    // storing information about input parameters
    FILE *input;
    input = filesaving::open_file("w+", params->output_directory_name, file.c_str());
    fprintf(input, "\nRelativistic Vlasov-Maxwell solver. Developed by Artem Korzhimanov. mailto: korzhimanov.artem@gmail.com\n");

    fprintf(input, "\nGENERAL PARAMETERS\n\n");
    fprintf(input, "\tSpace cells per wavelength = %d\n", params->mesh->ppw);
    fprintf(input, "\tSpace step = %f Wavelengths\n", .5*M_1_PI*params->mesh->dz);
    fprintf(input, "\tTotal space cells = %d\n", params->mesh->MAX_Z);
    fprintf(input, "\tTotal space size = %f Wavelengths\n", (double)params->mesh->MAX_Z/params->mesh->ppw);
    fprintf(input, "\tTime step = %f Waveperiods\n", .5*M_1_PI*params->mesh->dt);
    fprintf(input, "\tTotal time steps = %d\n", params->mesh->MAX_T);
    fprintf(input, "\tTotal running time = %f Waveperiods\n", .5*M_1_PI*params->mesh->MAX_T*params->mesh->dt);

    fprintf(input, "\nPLASMA PARAMETERS\n\n");
    for (int sp = 0; sp < plasmas->species_number; sp++)
    {
        fprintf(input, "\n\tSpecies %d\n\n", sp);
        plasmas->pfc[sp].SaveInput(input);
    }

    fprintf(input, "\nTEST PARTICLES PARAMETERS\n\n");
    fprintf(input, "\tNumber of particles = %d\n", particles->particles_number);
    if (particles->particles_number > 0)
    {
        fprintf(input, "\tMass of particles = %f\n", particles->mass);
        fprintf(input, "\tCharge of particles = %f\n", particles->charge);
        fprintf(input, "\tParticles are equidistantly deposited from %f to %f Wavelenths\n", (float)particles->start_point/params->mesh->ppw, (float)(particles->start_point + particles->interval*particles->particles_number)/params->mesh->ppw);
    }

    fprintf(input, "\nOUTPUT PARAMETERS\n\n");
    if (params->save_fields == true)
        fprintf(input, "\tFields are being saved every %d steps or every %f Waveperiods in files in %s-format\n", params->save_fields_dt, (double)params->save_fields_dt/params->mesh->ppw, params->save_fields_format.c_str());
    else fprintf(input, "\tFields are not being saved\n");
    if (params->save_concs == true)
        fprintf(input, "\tConcentration distributions are being saved every %d steps or every %f Waveperiods in files in %s-format\n", params->save_concs_dt, (double)params->save_concs_dt/params->mesh->ppw, params->save_concs_format.c_str());
    else fprintf(input, "\tConcentration distributions are not being saved\n");
    if (params->save_dstr == true)
        fprintf(input, "\tDistribution functions are being saved every %d steps or every %f Waveperiods in files in %s-format\n", params->save_dstr_dt, (double)params->save_dstr_dt/params->mesh->ppw, params->save_dstr_format.c_str());
    else fprintf(input, "\tDistribution functions are not being saved\n");

    filesaving::close_file(input);
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
        z = int(floor(particles->prt[i].r[2]/params->mesh->dz));
        if (z < 0 || z >= params->mesh->MAX_Z) continue;
        dz = particles->prt[i].r[2]/params->mesh->dz - floor(particles->prt[i].r[2]/params->mesh->dz);

        E[0] = (1. - dz) * fdtd->ex[z] + dz * fdtd->ex[z+1];
        E[1] = (1. - dz) * fdtd->ey[z] + dz * fdtd->ey[z+1];

        z = int(floor(particles->prt[i].r[2]/params->mesh->dz - 0.5));
        if (z < 0 || z >= params->mesh->MAX_Z-1) continue;
        dz = particles->prt[i].r[2]/params->mesh->dz - floor(particles->prt[i].r[2]/params->mesh->dz);

        E[2] = (1. - dz) * plasmas->ez[z] + dz * plasmas->ez[z+1];
        B[0] = (1. - dz) * fdtd->hx[z] + dz * fdtd->hx[z+1];
        B[1] = (1. - dz) * fdtd->hy[z] + dz * fdtd->hy[z+1];

        particles->prt[i].MakeStep(params->mesh->dt, E, B);
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
            filesaving::save_file_1D(plasmas->ez, params->mesh->MAX_Z, params->output_directory_name, "ez.txt");

            filesaving::save_file_1D(plasmas->a2, params->mesh->MAX_Z, params->output_directory_name, "vecpot2.txt");
        }

        if (format == "bin")
        {
            filesaving::save_file_1D_bin(fdtd->ex, params->mesh->MAX_Z, params->output_directory_name, "ex.bin");
            filesaving::save_file_1D_bin(fdtd->ey, params->mesh->MAX_Z, params->output_directory_name, "ey.bin");
            filesaving::save_file_1D_bin(fdtd->hx, params->mesh->MAX_Z, params->output_directory_name, "hx.bin");
            filesaving::save_file_1D_bin(fdtd->hy, params->mesh->MAX_Z, params->output_directory_name, "hy.bin");
            filesaving::save_file_1D_bin(plasmas->ez, params->mesh->MAX_Z, params->output_directory_name, "ez.bin");

            filesaving::save_file_1D_bin(plasmas->a2, params->mesh->MAX_Z, params->output_directory_name, "vecpot2.bin");
        }

        if (format == "gzip")
        {
            filesaving::save_file_1D_gzip(fdtd->ex, params->mesh->MAX_Z, params->output_directory_name, "ex.gz");
            filesaving::save_file_1D_gzip(fdtd->ey, params->mesh->MAX_Z, params->output_directory_name, "ey.gz");
            filesaving::save_file_1D_gzip(fdtd->hx, params->mesh->MAX_Z, params->output_directory_name, "hx.gz");
            filesaving::save_file_1D_gzip(fdtd->hy, params->mesh->MAX_Z, params->output_directory_name, "hy.gz");
            filesaving::save_file_1D_gzip(plasmas->ez, params->mesh->MAX_Z, params->output_directory_name, "ez.gz");

            filesaving::save_file_1D_gzip(plasmas->a2, params->mesh->MAX_Z, params->output_directory_name, "vecpot2.gz");
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

        for (int sp = 0; sp < plasmas->species_number; sp++)
        {
            std::stringstream ss;
            ss << params->output_directory_name << "conc" << sp;
            if (format == "txt")
                plasmas->pfc[sp].SaveConcentrationTxt(ss.str());

            if (format == "bin")
                plasmas->pfc[sp].SaveConcentrationBin(ss.str());

            if (format == "gzip")
                plasmas->pfc[sp].SaveConcentrationGZip(ss.str());
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

        for (int sp = 0; sp < plasmas->species_number; sp++)
        {
            std::stringstream ss;
            ss << params->output_directory_name << "dstr_func" << sp << "/data" << std::setfill('0') << std::setw(8) << k;
            if (format == "txt")
                plasmas->pfc[sp].SaveDstrFunctionTxt(ss.str());

            if (format == "bin")
                plasmas->pfc[sp].SaveDstrFunctionBin(ss.str());

            if (format == "gzip")
                plasmas->pfc[sp].SaveDstrFunctionGZip(ss.str());
        }
    }
}

void Solver::SavePrtDat(int k)
{
    if (particles->particles_number > 0)
    {
        std::ofstream *out_file;
        for (int i = 0; i < particles->particles_number; i++)
        {
            char file_path[256];
            strcpy(file_path, params->output_directory_name.c_str());
            char file_name[256];
            sprintf(file_name, "particles/particle%06d.txt", i);
            strcat(file_path, file_name);
            out_file = new std::ofstream(file_path, std::ios_base::app);
            *out_file << particles->prt[i].r[0] << "\t" << particles->prt[i].r[1] << "\t" << particles->prt[i].r[2] << "\t" << particles->prt[i].p[0] << "\t" << particles->prt[i].p[1] << "\t" << particles->prt[i].p[2] << std::endl;
            out_file->close();
            delete out_file;
        }
    }
}

void Solver::SaveOutput(int k, std::string file)
{
    if (k%params->mesh->ppw == 0)
    {
        for (int sp = 0; sp < plasmas->species_number; sp++)
            plasma_energy[sp] = plasmas->pfc[sp].KineticEnergy(plasmas->ax, plasmas->ay);

        laser_energy = fdtd->Energy();

        electrostatic_energy = 0.;
        for (int i = 0; i <= params->mesh->MAX_Z; i++)
            electrostatic_energy += plasmas->ez[i]*plasmas->ez[i];
        electrostatic_energy *= 0.5*params->mesh->dz;

        FILE *output;
        output = filesaving::open_file("a+", params->output_directory_name, file.c_str());
        fprintf(output, "%06d|%16.8e|", k, plasmas->ez[params->mesh->MAX_Z-1]);
        for (int sp = 0; sp < plasmas->species_number; sp++)
            fprintf(output, " %16.8e|", plasma_energy[sp]);
        fprintf(output, "%16.8e|%16.8e|%16.8e|%16.8e|\n", electrostatic_energy, laser_energy, em_flux, mymath::sum(plasma_energy, plasmas->species_number)+electrostatic_energy+laser_energy-em_flux);
        filesaving::close_file(output);

        if (k==0) plasma_energy_0 = mymath::sum(plasma_energy, plasmas->species_number);
    }
}

void Solver::SaveResults()
{
    FILE *en_out;
    en_out = filesaving::open_file("a+", "", "energy.txt");
    double plasma_energy_1 = mymath::sum(plasma_energy, plasmas->species_number)+electrostatic_energy+laser_energy;
    fprintf(en_out, "%.8g\t%.8g\t%.8g\n", plasma_energy_0, plasma_energy_1, 2.*(plasma_energy_1-plasma_energy_0));
    filesaving::close_file(en_out);

    en_out = filesaving::open_file("a+", "", "energy2.txt");
    fprintf(en_out, "%.8g\t", 2.*(plasma_energy_1-plasma_energy_0));
    filesaving::close_file(en_out);
}

#endif // SOLVER_H
