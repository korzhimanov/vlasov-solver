/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file output.cpp
 * \brief The file implements methods of Output class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 *
 * The file contains argument parser, timer, main cycle and basic logging.
 */


#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include "output.h"
#include "file_saving.h"
#include "solver.h"


Output::Output(pyinput* in, std::string dn, Solver* s, int* err) : save_fields(1),
                                                         save_concs(1),
                                                         save_dstr(0),
                                                         save_dt(1),
                                                         save_fields_dt(0),
                                                         save_concs_dt(0),
                                                         save_dstr_dt(0),
                                                         save_format("gzip"),
                                                         save_fields_format(""),
                                                         save_concs_format(""),
                                                         save_dstr_format(""),
                                                         solver(s)
{
    *err = Init(in, &dn);
    CreateDirs();
    SaveInput("initial_parameters.txt");
    InitEnergyFile("energy.txt");
}

Output::~Output()
{
    if (plasma_energy) delete[] plasma_energy;
    if (energy_file) fclose(energy_file);
}

int Output::Init(pyinput* in, std::string *dn)
{
    output_directory_name = *dn + "/";

    save_fields = in->GetInt("save_fields");
    save_concs  = in->GetInt("save_concs");
    save_dstr   = in->GetInt("save_dstr");

    if ( (save_fields || save_concs || save_dstr) == true)
    {
        if ( !in->SetPositive("save_dt", &save_dt) ) return 1010;
        save_format = "txt"; if ( !SetFormat(in, "save_format", &(save_format)) ) return 1020;

        if (save_fields == true)
        {
            if ( !in->SetNotNegative("save_fields_dt", &save_fields_dt) ) return 1030;
            if (save_fields_dt == 0)
                save_fields_dt = save_dt;
            if ( !SetFormat(in, "save_fields_format", &(save_fields_format)) ) return 1031;
            if (save_fields_format == "")
                save_fields_format = save_format;
        }
        if (save_concs == true)
        {
            if ( !in->SetNotNegative("save_concs_dt", &save_concs_dt) ) return 1040;
            if (save_concs_dt == 0)
                save_concs_dt = save_dt;
            if ( !SetFormat(in, "save_concs_format", &(save_concs_format)) ) return 1041;
            if (save_concs_format == "")
                save_concs_format = save_format;
        }
        if (save_dstr == true)
        {
            if ( !in->SetNotNegative("save_dstr_dt", &save_dstr_dt) ) return 1050;
            if (save_dstr_dt == 0)
                save_dstr_dt = save_dt;
            if ( !SetFormat(in, "save_dstr_format", &(save_dstr_format)) ) return 1051;
            if (save_dstr_format == "")
                save_dstr_format = save_format;
        }
    }

    return 0;
}

bool Output::SetFormat(pyinput* in, std::string name, std::string *var)
{
    *var = in->GetString(name);
    if (save_format == "txt" || save_format == "bin" || save_format ==  "gzip") return true;
    else {std::cout << name + " must be either 'txt', 'bin' or 'gzip'" << std::endl; return false;}
}

void Output::CreateDirs()
{
    // creating main directory for output
    {
        filesaving::create_dir("", (output_directory_name).c_str());
    }

    // creating directory for test particles
    if (solver->particles->particles_number > 0)
    {
        filesaving::create_dir(output_directory_name, "particles/");
    }

    // creating directories for storing distribution functions in files
    if (save_dstr == true)
        for (int sp = 0; sp < solver->plasmas->species_number; sp++)
        {
            filesaving::create_dir(output_directory_name, "dstr_func%d/", sp);
        }
}

void Output::SaveInput(std::string fn)
{
    // storing information about input parameters
    FILE *input;
    input = fopen((output_directory_name + fn).c_str(), "w+");
    fprintf(input, "\nRelativistic Vlasov-Maxwell solver. Developed by Artem Korzhimanov. mailto: korzhimanov.artem@gmail.com\n");

    fprintf(input, "\nGENERAL PARAMETERS\n\n");
    fprintf(input, "\tSpace cells per wavelength = %d\n", solver->mesh->ppw);
    fprintf(input, "\tSpace step = %f Wavelengths\n", .5*M_1_PI*solver->mesh->dz);
    fprintf(input, "\tTotal space cells = %d\n", solver->mesh->MAX_Z);
    fprintf(input, "\tTotal space size = %f Wavelengths\n", (double)solver->mesh->MAX_Z/solver->mesh->ppw);
    fprintf(input, "\tTime step = %f Waveperiods\n", .5*M_1_PI*solver->mesh->dt);
    fprintf(input, "\tTotal time steps = %d\n", solver->mesh->MAX_T);
    fprintf(input, "\tTotal running time = %f Waveperiods\n", .5*M_1_PI*solver->mesh->MAX_T*solver->mesh->dt);

    fprintf(input, "\nPLASMA PARAMETERS\n\n");
    for (int sp = 0; sp < solver->plasmas->species_number; sp++)
    {
        fprintf(input, "\n\tSpecies %d\n\n", sp);
        solver->plasmas->pfc[sp].SaveInput(input);
    }

    fprintf(input, "\nTEST PARTICLES PARAMETERS\n\n");
    fprintf(input, "\tNumber of particles = %d\n", solver->particles->particles_number);
    if (solver->particles->particles_number > 0)
    {
        fprintf(input, "\tMass of particles = %f\n", solver->particles->mass);
        fprintf(input, "\tCharge of particles = %f\n", solver->particles->charge);
        fprintf(input, "\tParticles are equidistantly deposited from %f to %f Wavelenths\n", (float)solver->particles->start_point/solver->mesh->ppw, (float)(solver->particles->start_point + solver->particles->interval*solver->particles->particles_number)/solver->mesh->ppw);
    }

    fprintf(input, "\nOUTPUT PARAMETERS\n\n");
    if (save_fields == true)
        fprintf(input, "\tFields are being saved every %d steps or every %f Waveperiods in files in %s-format\n", save_fields_dt, (double)save_fields_dt/solver->mesh->ppw, save_fields_format.c_str());
    else fprintf(input, "\tFields are not being saved\n");
    if (save_concs == true)
        fprintf(input, "\tConcentration distributions are being saved every %d steps or every %f Waveperiods in files in %s-format\n", save_concs_dt, (double)save_concs_dt/solver->mesh->ppw, save_concs_format.c_str());
    else fprintf(input, "\tConcentration distributions are not being saved\n");
    if (save_dstr == true)
        fprintf(input, "\tDistribution functions are being saved every %d steps or every %f Waveperiods in files in %s-format\n", save_dstr_dt, (double)save_dstr_dt/solver->mesh->ppw, save_dstr_format.c_str());
    else fprintf(input, "\tDistribution functions are not being saved\n");

    fclose(input);
}

void Output::SaveFields(int k)
{
    em_flux += solver->fdtd->FluxIn();

    if (save_fields == true && k%save_fields_dt == 0)
    {
        if (save_fields_format == "txt")
        {
            filesaving::save_file_1D(solver->fdtd->ex, solver->mesh->MAX_Z, output_directory_name, "ex.txt");
            filesaving::save_file_1D(solver->fdtd->ey, solver->mesh->MAX_Z, output_directory_name, "ey.txt");
            filesaving::save_file_1D(solver->fdtd->hx, solver->mesh->MAX_Z, output_directory_name, "hx.txt");
            filesaving::save_file_1D(solver->fdtd->hy, solver->mesh->MAX_Z, output_directory_name, "hy.txt");
            filesaving::save_file_1D(solver->plasmas->ez, solver->mesh->MAX_Z, output_directory_name, "ez.txt");

            filesaving::save_file_1D(solver->plasmas->a2, solver->mesh->MAX_Z, output_directory_name, "vecpot2.txt");
        }

        if (save_fields_format == "bin")
        {
            filesaving::save_file_1D_bin(solver->fdtd->ex, solver->mesh->MAX_Z, output_directory_name, "ex.bin");
            filesaving::save_file_1D_bin(solver->fdtd->ey, solver->mesh->MAX_Z, output_directory_name, "ey.bin");
            filesaving::save_file_1D_bin(solver->fdtd->hx, solver->mesh->MAX_Z, output_directory_name, "hx.bin");
            filesaving::save_file_1D_bin(solver->fdtd->hy, solver->mesh->MAX_Z, output_directory_name, "hy.bin");
            filesaving::save_file_1D_bin(solver->plasmas->ez, solver->mesh->MAX_Z, output_directory_name, "ez.bin");

            filesaving::save_file_1D_bin(solver->plasmas->a2, solver->mesh->MAX_Z, output_directory_name, "vecpot2.bin");
        }

        if (save_fields_format == "gzip")
        {
            filesaving::save_file_1D_gzip(solver->fdtd->ex, solver->mesh->MAX_Z, output_directory_name, "ex.gz");
            filesaving::save_file_1D_gzip(solver->fdtd->ey, solver->mesh->MAX_Z, output_directory_name, "ey.gz");
            filesaving::save_file_1D_gzip(solver->fdtd->hx, solver->mesh->MAX_Z, output_directory_name, "hx.gz");
            filesaving::save_file_1D_gzip(solver->fdtd->hy, solver->mesh->MAX_Z, output_directory_name, "hy.gz");
            filesaving::save_file_1D_gzip(solver->plasmas->ez, solver->mesh->MAX_Z, output_directory_name, "ez.gz");

            filesaving::save_file_1D_gzip(solver->plasmas->a2, solver->mesh->MAX_Z, output_directory_name, "vecpot2.gz");
        }

    }
}

void Output::SaveConcs(int k)
{
    if (save_concs == true && k%save_concs_dt == 0)
    {
        for (int sp = 0; sp < solver->plasmas->species_number; sp++)
        {
            std::stringstream ss;
            ss << output_directory_name << "conc" << sp;
            if (save_concs_format == "txt")
                solver->plasmas->pfc[sp].SaveConcentrationTxt(ss.str());

            if (save_concs_format == "bin")
                solver->plasmas->pfc[sp].SaveConcentrationBin(ss.str());

            if (save_concs_format == "gzip")
                solver->plasmas->pfc[sp].SaveConcentrationGZip(ss.str());
        }
    }
}

void Output::SaveDstrFunc(int k)
{
    if (save_dstr == true && k%save_dstr_dt == 0)
    {
        for (int sp = 0; sp < solver->plasmas->species_number; sp++)
        {
            std::stringstream ss;
            ss << output_directory_name << "dstr_func" << sp << "/data" << std::setfill('0') << std::setw(8) << k;
            if (save_dstr_format == "txt")
                solver->plasmas->pfc[sp].SaveDstrFunctionTxt(ss.str());

            if (save_dstr_format == "bin")
                solver->plasmas->pfc[sp].SaveDstrFunctionBin(ss.str());

            if (save_dstr_format == "gzip")
                solver->plasmas->pfc[sp].SaveDstrFunctionGZip(ss.str());
        }
    }
}

void Output::SavePrtData(int k)
{
    if (solver->particles->particles_number > 0)
    {
        std::ofstream *out_file;
        for (int i = 0; i < solver->particles->particles_number; i++)
        {
            char file_path[256];
            strcpy(file_path, output_directory_name.c_str());
            char file_name[256];
            sprintf(file_name, "particles/particle%06d.txt", i);
            strcat(file_path, file_name);
            out_file = new std::ofstream(file_path, std::ios_base::app);
            *out_file << solver->particles->prt[i].r[0] << "\t" << solver->particles->prt[i].r[1] << "\t" << solver->particles->prt[i].r[2] << "\t" << solver->particles->prt[i].p[0] << "\t" << solver->particles->prt[i].p[1] << "\t" << solver->particles->prt[i].p[2] << std::endl;
            out_file->close();
            delete out_file;
        }
    }
}

void Output::InitEnergyFile(std::string fn)
{
    em_flux = 0.;
    plasma_energy = new double[solver->plasmas->species_number];
    energy_file = fopen((output_directory_name + fn).c_str(), "w+");
    fprintf(energy_file, " Step |   Neutrality   |");
    for (int sp = 0; sp < solver->plasmas->species_number; sp++)
        fprintf(energy_file, "       W_%d       |", sp);
    fprintf(energy_file, "     W_el-st    |     W_laser    |     EM Flux    |     W_total    |     Balance    |\n");
    fflush(energy_file);
}

void Output::WriteEnergy(int k)
{
    if (k%solver->mesh->ppw == 0)
    {
        for (int sp = 0; sp < solver->plasmas->species_number; sp++)
            plasma_energy[sp] = solver->plasmas->pfc[sp].KineticEnergy(solver->plasmas->ax, solver->plasmas->ay);

        if (k==0) plasma_energy_0 = mymath::sum(plasma_energy, solver->plasmas->species_number);

        laser_energy = solver->fdtd->Energy();

        electrostatic_energy = 0.;
        for (int i = 0; i <= solver->mesh->MAX_Z; i++)
            electrostatic_energy += solver->plasmas->ez[i]*solver->plasmas->ez[i];
        electrostatic_energy *= 0.5*solver->mesh->dz;

        fprintf(energy_file, "%06d|%16.8e|", k, solver->plasmas->ez[solver->mesh->MAX_Z-1]);
        for (int sp = 0; sp < solver->plasmas->species_number; sp++)
            fprintf(energy_file, " %16.8e|", plasma_energy[sp]);

        total_energy = mymath::sum(plasma_energy, solver->plasmas->species_number)+electrostatic_energy+laser_energy;
        fprintf(energy_file, "%16.8e|%16.8e|%16.8e|%16.8e|%16.8e|\n", electrostatic_energy, laser_energy, em_flux, total_energy, total_energy-em_flux-plasma_energy_0);
        fflush(energy_file);
    }
}
