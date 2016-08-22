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
#include "pyinput.h"
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
        FileSaving *fs;
        MyMath *mm;
        // general parameters
        int ppw; Mesh *mesh;
        double THETA;
        // plasma parameters
        int NUM_SP;
        double N_0;
        // pulse parameters
        pFunc *pulse_x, *pulse_y;
        // output parameters
        string output_directory_name;
        bool save_fields, save_concs, save_dstr;
        int save_dt, save_fields_dt, save_concs_dt, save_dstr_dt;
        string save_format, save_fields_format, save_concs_format, save_dstr_format;
        // particles
        Particle *prt; int NUM_PRT, start_point; double MASS_PRT, CHARGE_PRT, interval;

        double *fixed_ions_conc; // concentration of fixed ions
        pFunc *fixed_ions_profile; // fixed ions concentration profile function
        FDTD *fdtd; // FDTD - class
        PFC *pfc; // PFC - class
        int source; // position of a source of an electromagnetic wave in respect to a PML-layer
        double *ax, *ay, *a2; // vector potential
        double *ez; // longitudinal field

        char *main_dir;
        FILE *output;

        pyinput *in;

        double em_flux; // flux density of electromagnetic field
        double electrostatic_energy, laser_energy, *plasma_energy, plasma_energy_0;

    public:
//------constructors---------------------------------------------------
        Solver();
        Solver(string file);
//------destructor-----------------------------------------------------
        ~Solver();

//------initialisation-------------------------------------------------
        int InitVars(string file);
        void AllocMemory();
        void InitFields();
        void InitPlasma(); // initializes plasma distribution
        void InitTestParticles();
        void CreateDirs();
        void SaveInput(string file);
        void InitOutput(string file);

//------calculations---------------------------------------------------
        void MoveParticles();
        inline void CalcSources(double t)
        {
            fdtd->exl =   pulse_x->call((double)(t - .5*mesh->dz));
            fdtd->eyl =   pulse_y->call((double)(t - .5*mesh->dz));
            fdtd->hxl = - pulse_y->call((double)(t + .5*mesh->dt));
            fdtd->hyl =   pulse_x->call((double)(t + .5*mesh->dt));
        }

        inline void CalcLongFields()
        {
            // calculates longitudinal field
            memset(ez, 0, sizeof(double)*(mesh->MAX_Z+1));
            for (int sp = 0; sp < NUM_SP; sp++)
                pfc[sp].CalcLongitudinalField(ez);

            ez[0] += 0.5*N_0*mesh->dz*fixed_ions_conc[0];
            for (int i = 1; i < mesh->MAX_Z; i++)
                ez[i] += ez[i-1] + 0.5*N_0*mesh->dz*(fixed_ions_conc[i]+fixed_ions_conc[i-1]);
            ez[mesh->MAX_Z] += 0.5*N_0*mesh->dz*fixed_ions_conc[mesh->MAX_Z-1];
        }

        inline void CalcTransFields()
        {
            // evaluates vector potential
            for (int i = 0; i <= mesh->MAX_Z; i++)
            {
                ax[i] -= fdtd->ex[i]*mesh->dt;
                ay[i] -= fdtd->ey[i]*mesh->dt;
                a2[i] = ax[i]*ax[i]+ay[i]*ay[i];
            }

            // evaluating electric and magnetic fields
            fdtd->Maxwell();

            em_flux += fdtd->FluxIn();
        }

        inline void FieldGeneration()
        {
            // field generation by plasma currents
            for (int sp = 0; sp < NUM_SP; sp++)
                pfc[sp].CalcCurrent(fdtd, ax, ay, a2);

            fdtd->ey[0] -= 0.5*sin(THETA)*N_0*mesh->dz*fixed_ions_conc[0];
            for (int i = 1; i < mesh->MAX_Z; i++)
                fdtd->ey[i] -= 0.5*sin(THETA)*N_0*mesh->dz*(fixed_ions_conc[i]+fixed_ions_conc[i-1]);
            fdtd->ey[mesh->MAX_Z] -= 0.5*sin(THETA)*N_0*mesh->dz*fixed_ions_conc[mesh->MAX_Z-1];
        }

        inline void CalcDstrFunc()
        {
            for (int sp = 0; sp < NUM_SP; sp++)
                pfc[sp].MakeStep(ez, ax, ay);
        }

//------saving data----------------------------------------------------
        void SaveFields(int k); // saving fields data in files
        void SaveConcs(int k); // saving concentarations data in files
        void SaveDstrFunc(int k); // saving distribution functions data in files
        void SavePrtDat(int k); // save particles data
        void SaveOutput(int k, string file);
        void SaveResults();

    private:
//------miscellaneous--------------------------------------------------
        bool SetNotNegative(pyinput *in, string name, int *var);
        bool SetPositive(pyinput *in, string name, int *var);
        bool SetPositive(pyinput *in, string name, double *var);
        bool SetFormat(pyinput *in, string name, string *var);
};

Solver::Solver(string file = "input.py")
{
    fs = new FileSaving;
    mm = new MyMath;

    int err = InitVars(file);
    if( err == 0 )
    {
        cout << "Initializing succesful!" << endl;
        AllocMemory();
    }
    else
        cout << "Error! Initializing failed! Returning code is " << err << endl;
}

Solver::~Solver()
{
    delete fixed_ions_profile;
    delete pulse_x;
    delete pulse_y;

    delete mesh;

    delete[] pfc;

    delete[] fixed_ions_conc;

    delete fdtd;
    delete[] ax;
    delete[] ay;
    delete[] a2;
    delete[] ez;

    if (NUM_PRT > 0) delete[] prt;

    delete[] main_dir;
    delete[] plasma_energy;

    delete in;

    delete mm;
    delete fs;
}

int Solver::InitVars(string file)
{
    mesh = new Mesh;
    in = new pyinput;

    cout << "Reading input file " << file << " ... " << endl;
    in->ReadFile(file);
    cout << "done!" << endl;

    int tpi;
    double tpf;

    ppw    = 16; if ( SetPositive(in,   "ppw",    &ppw) ) return 100;
    tpi   = 32; if ( SetPositive(in, "MAX_Z",     &tpi) ) return 200; mesh->MAX_Z = tpi;
    tpi   = 32; if ( SetNotNegative(in, "MAX_T",  &tpi) ) return 210; mesh->MAX_T = tpi;
    NUM_SP =  1; if ( SetPositive(in,"NUM_SP", &NUM_SP) ) return 300;
    tpf = 2.*M_PI/ppw; if ( SetPositive(in, "dz", &tpf) ) return 400; mesh->dz = tpf;
    tpf = 2.*M_PI/ppw; if ( SetPositive(in, "dt", &tpf) ) return 410; mesh->dt = tpf;
    mesh->dt_dz = mesh->dt/mesh->dz;
    THETA = in->GetDouble("THETA");

    pfc = new PFC[NUM_SP];
    for (int sp = 0; sp < NUM_SP; sp++)
        pfc[sp].Init(sp, in, mesh);

    fixed_ions_profile = new pFunc(in->GetFunc("FIXED_IONS_PROFILE"));

    N_0 = 0.; if ( SetPositive(in,"N_0", &N_0) ) return 1200;
    if (THETA >= 0. && THETA < M_PI/2)
        N_0 = N_0/(cos(THETA)*cos(THETA)*cos(THETA));
    else
        cout << "Invalid incident angle! It should be between 0 and pi/2." << endl;

    pulse_x = new pFunc(in->GetFunc("PULSE_X"));
    pulse_y = new pFunc(in->GetFunc("PULSE_Y"));
    source = 1; if ( SetPositive(in, "source", &source) ) return 1400;

    NUM_PRT = 0; if ( SetNotNegative(in, "NUM_PRT", &NUM_PRT) ) return  1600;
    if (NUM_PRT > 0)
    {
        start_point =  0; if ( SetNotNegative(in, "start_point", &start_point) ) return  1610;
        interval    = 1.; if ( SetPositive(in,    "interval", &interval   ) ) return  1620;
        MASS_PRT    = 1.; if ( SetPositive(in,    "MASS_PRT", &MASS_PRT   ) ) return  1630;
        CHARGE_PRT  = in->GetDouble("CHARGE_PRT");
    }

    output_directory_name = in->GetString("output_directory_name");

    save_fields = in->GetInt("save_fields");
    save_concs  = in->GetInt("save_concs");
    save_dstr   = in->GetInt("save_dstr");

    if ( (save_fields || save_concs || save_dstr) == true)
    {
        save_dt     =   ppw; if ( SetPositive(in,     "save_dt", &save_dt)    ) return 1410;
        save_format = "txt"; if (   SetFormat(in, "save_format", &save_format)) return 1420;

        if (save_fields == true)
        {
            save_fields_dt     =  0; if ( SetNotNegative(in,     "save_fields_dt", &save_fields_dt)    ) return 1430;
            if (save_fields_dt == 0)
                save_fields_dt = save_dt;
            save_fields_format = ""; if (      SetFormat(in, "save_fields_format", &save_fields_format)) return 1431;
            if (save_fields_format == "")
                save_fields_format = save_format;
        }
        if (save_concs == true)
        {
            save_concs_dt     =  0; if ( SetNotNegative(in,     "save_concs_dt", &save_concs_dt)    ) return 1440;
            if (save_concs_dt == 0)
                save_concs_dt = save_dt;
            save_concs_format = ""; if (      SetFormat(in, "save_concs_format", &save_concs_format)) return 1441;
            if (save_concs_format == "")
                save_concs_format = save_format;
        }
        if (save_dstr == true)
        {
            save_dstr_dt     =  0; if ( SetNotNegative(in,     "save_dstr_dt", &save_dstr_dt)    ) return 1450;
            if (save_dstr_dt == 0)
                save_dstr_dt = save_dt;
            save_dstr_format = ""; if (      SetFormat(in, "save_dstr_format", &save_dstr_format)) return 1451;
            if (save_dstr_format == "")
                save_dstr_format = save_format;
        }
    }

    return 0;
}

void Solver::AllocMemory()
{
    // memory allocation for fixed ions
    fixed_ions_conc = new double[mesh->MAX_Z];
    for (int sp = 0; sp < NUM_SP; sp++)
        pfc[sp].AllocMemory();

    // memory allocation for electromagnetic fields
    fdtd = new FDTD;
    ax = new double[mesh->MAX_Z + 1];
    ay = new double[mesh->MAX_Z + 1];
    a2 = new double[mesh->MAX_Z + 1];
    ez = new double[mesh->MAX_Z + 1];

    // memory allocation for test particles
    if (NUM_PRT > 0) prt = new Particle[NUM_PRT];

    main_dir = new char[512];
    plasma_energy = new double[NUM_SP];
}

void Solver::InitFields()
{
    fdtd->Init(mesh, source);
    fdtd->InitPML(int(1./mesh->dz), 10.);
    mm->zeros(ax, mesh->MAX_Z+1);
    mm->zeros(ay, mesh->MAX_Z+1);
    mm->zeros(a2, mesh->MAX_Z+1);
    mm->zeros(ez, mesh->MAX_Z+1);

    em_flux = 0.;
}

void Solver::InitPlasma()
{
    for (int i = 0; i < mesh->MAX_Z; i++)
        fixed_ions_conc[i] = fixed_ions_profile->call(i);
    for (int sp = 0; sp < NUM_SP; sp++)
        pfc[sp].SetDistribution();
}

void Solver::InitTestParticles()
{
    if (NUM_PRT > 0) for (int i = 0; i < NUM_PRT; i++)
    {
        prt[i].q = int(CHARGE_PRT);
        prt[i].m = MASS_PRT;
        for (int j = 0; j < 3; j++)
        {
            prt[i].r[j] = 0.;
            prt[i].p[j] = 0.;
        }
        prt[i].p[1] = tan(THETA);
        prt[i].r[2] = (start_point + interval*i)*mesh->dz;
    }
}

void Solver::CreateDirs()
{
    // creating main directory for saving
    if (output_directory_name == "")
    {
        if (SYSTEM == 0)
            fs->create_dir(main_dir, "", "n0_%.3g\\", N_0);
        else
            fs->create_dir(main_dir, "", "n0_%.3g/", N_0);
    }
    else
    {
        if (SYSTEM == 0)
            fs->create_dir(main_dir, "", (output_directory_name + "\\").c_str());
        else
            fs->create_dir(main_dir, "", (output_directory_name + "/").c_str());
    }

    // creating directory for test particles
    if (NUM_PRT > 0)
    {
        if (SYSTEM == 0)
            fs->create_dir(main_dir, "particles\\");
        else
            fs->create_dir(main_dir, "particles/");
    }

    // creating directories for storing distribution functions in files
    if (save_dstr == true)
        for (int sp = 0; sp < NUM_SP; sp++)
        {
            if (SYSTEM == 0)
                fs->create_dir(main_dir, "dstr_func%d\\", sp);
            else
                fs->create_dir(main_dir, "dstr_func%d/", sp);
        }
}

void Solver::SaveInput(string file)
{
    // storing information about input parameters
    FILE *input;
    input = fs->open_file("w+", main_dir, file.c_str());
    fprintf(input, "\nVlasov-Maxwell code. Developed by Artem Korzhimanov and Arkady Gonoskov. mailto: kav@ufp.appl.sci-nnov.ru\n");

    fprintf(input, "\nGENERAL PARAMETERS\n\n");
    fprintf(input, "\tSpace cells per wavelength = %d\n", ppw);
    fprintf(input, "\tSpace step = %f Wavelengths\n", .5*M_1_PI*mesh->dz);
    fprintf(input, "\tTotal space cells = %d\n", mesh->MAX_Z);
    fprintf(input, "\tTotal space size = %f Wavelengths\n", (double)mesh->MAX_Z/ppw);
    fprintf(input, "\tTime step = %f Waveperiods\n", .5*M_1_PI*mesh->dt);
    fprintf(input, "\tTotal time steps = %d\n", mesh->MAX_T);
    fprintf(input, "\tTotal running time = %f Waveperiods\n", .5*M_1_PI*mesh->MAX_T*mesh->dt);

    fprintf(input, "\nPLASMA PARAMETERS\n\n");
    fprintf(input, "\tOverdense parameter = %f\n", N_0);
    for (int sp = 0; sp < NUM_SP; sp++)
    {
        fprintf(input, "\n\tSpecies %d\n\n", sp);
        pfc[sp].SaveInput(input);
    }

    fprintf(input, "\nTEST PARTICLES PARAMETERS\n\n");
    fprintf(input, "\tNumber of particles = %d\n", NUM_PRT);
    if (NUM_PRT > 0)
    {
        fprintf(input, "\tMass of particles = %f\n", MASS_PRT);
        fprintf(input, "\tCharge of particles = %f\n", CHARGE_PRT);
        fprintf(input, "\tParticles are equidistantly deposited from %f to %f Wavelenths\n", (float)start_point/ppw, (float)(start_point + interval*NUM_PRT)/ppw);
    }

    fprintf(input, "\nOUTPUT PARAMETERS\n\n");
    if (save_fields == true)
        fprintf(input, "\tFields are being saved every %d steps or every %f Waveperiods in files in %s-format\n", save_fields_dt, (double)save_fields_dt/ppw, save_fields_format.c_str());
    else fprintf(input, "\tFields are not being saved\n");
    if (save_concs == true)
        fprintf(input, "\tConcentration distributions are being saved every %d steps or every %f Waveperiods in files in %s-format\n", save_concs_dt, (double)save_concs_dt/ppw, save_concs_format.c_str());
    else fprintf(input, "\tConcentration distributions are not being saved\n");
    if (save_dstr == true)
        fprintf(input, "\tDistribution functions are being saved every %d steps or every %f Waveperiods in files in %s-format\n", save_dstr_dt, (double)save_dstr_dt/ppw, save_dstr_format.c_str());
    else fprintf(input, "\tDistribution functions are not being saved\n");

    fs->close_file(input);
}

void Solver::InitOutput(string file)
{
    output = fs->open_file("w+", main_dir, file.c_str());
    fprintf(output, " Step |   Neutrality   |");
    for (int sp = 0; sp < NUM_SP; sp++)
        fprintf(output, "       W_%d       |", sp);
    fprintf(output, "     W_el-st    |     W_laser    |     EM Flux    |     W_total    |\n");
    fs->close_file(output);
}

void Solver::MoveParticles()
{
    int z;
    double dz;

    int i;

    double *E = new double[3];
    double *B = new double[3];
    B[2] = 0;
    for (i = 0; i < NUM_PRT; i++)
    {
        z = int(floor(prt[i].r[2]/mesh->dz));
        if (z < 0 || z >= mesh->MAX_Z) continue;
        dz = prt[i].r[2]/mesh->dz - floor(prt[i].r[2]/mesh->dz);

        E[0] = (1. - dz) * fdtd->ex[z] + dz * fdtd->ex[z+1];
        E[1] = (1. - dz) * fdtd->ey[z] + dz * fdtd->ey[z+1];

        z = int(floor(prt[i].r[2]/mesh->dz - 0.5));
        if (z < 0 || z >= mesh->MAX_Z-1) continue;
        dz = prt[i].r[2]/mesh->dz - floor(prt[i].r[2]/mesh->dz);

        E[2] = (1. - dz) * ez[z] + dz * ez[z+1];
        B[0] = (1. - dz) * fdtd->hx[z] + dz * fdtd->hx[z+1];
        B[1] = (1. - dz) * fdtd->hy[z] + dz * fdtd->hy[z+1];

        prt[i].MakeStep(mesh->dt, E, B);
    }

    delete[] E;
    delete[] B;
}

void Solver::SaveFields(int k)
{
    if (save_fields == true && k%save_fields_dt == 0)
    {
        string format;
        if (save_fields_format == "")
            format = save_format;
        else
            format = save_fields_format;

        if (format == "txt")
        {
            fs->save_file_1D(fdtd->ex, mesh->MAX_Z, main_dir, "ex.txt");
            fs->save_file_1D(fdtd->ey, mesh->MAX_Z, main_dir, "ey.txt");
            fs->save_file_1D(fdtd->hx, mesh->MAX_Z, main_dir, "hx.txt");
            fs->save_file_1D(fdtd->hy, mesh->MAX_Z, main_dir, "hy.txt");
            fs->save_file_1D(ez, mesh->MAX_Z, main_dir, "ez.txt");

            fs->save_file_1D(a2, mesh->MAX_Z, main_dir, "vecpot2.txt");
        }

        if (format == "bin")
        {
            fs->save_file_1D_bin(fdtd->ex, mesh->MAX_Z, main_dir, "ex.bin");
            fs->save_file_1D_bin(fdtd->ey, mesh->MAX_Z, main_dir, "ey.bin");
            fs->save_file_1D_bin(fdtd->hx, mesh->MAX_Z, main_dir, "hx.bin");
            fs->save_file_1D_bin(fdtd->hy, mesh->MAX_Z, main_dir, "hy.bin");
            fs->save_file_1D_bin(ez, mesh->MAX_Z, main_dir, "ez.bin");

            fs->save_file_1D_bin(a2, mesh->MAX_Z, main_dir, "vecpot2.bin");
        }

        if (format == "gzip")
        {
            fs->save_file_1D_gzip(fdtd->ex, mesh->MAX_Z, main_dir, "ex.gz");
            fs->save_file_1D_gzip(fdtd->ey, mesh->MAX_Z, main_dir, "ey.gz");
            fs->save_file_1D_gzip(fdtd->hx, mesh->MAX_Z, main_dir, "hx.gz");
            fs->save_file_1D_gzip(fdtd->hy, mesh->MAX_Z, main_dir, "hy.gz");
            fs->save_file_1D_gzip(ez, mesh->MAX_Z, main_dir, "ez.gz");

            fs->save_file_1D_gzip(a2, mesh->MAX_Z, main_dir, "vecpot2.gz");
        }

    }
}

void Solver::SaveConcs(int k)
{
    if (save_concs == true && k%save_concs_dt == 0)
    {
        string format;
        if (save_concs_format == "")
            format = save_format;
        else
            format = save_concs_format;

        for (int sp = 0; sp < NUM_SP; sp++)
        {
            stringstream ss;
            ss << main_dir << "conc" << sp;
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
    if (save_dstr == true && k%save_dstr_dt == 0)
    {
        string format;
        if (save_dstr_format == "")
            format = save_format;
        else
            format = save_dstr_format;

        for (int sp = 0; sp < NUM_SP; sp++)
        {
            stringstream ss;
            ss << main_dir << "dstr_func" << sp << "/data" << setfill('0') << setw(8) << k;
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
    if (NUM_PRT > 0)
    {
        ofstream *out_file;
        for (int i = 0; i < NUM_PRT; i++)
        {
            char file_path[256];
            strcpy(file_path, main_dir);
            char file_name[256];
            sprintf(file_name, "particles/particle%06d.txt", i);
            strcat(file_path, file_name);
            out_file = new ofstream(file_path, ios_base::app);
            *out_file << prt[i].r[0] << "\t" << prt[i].r[1] << "\t" << prt[i].r[2] << "\t" << prt[i].p[0] << "\t" << prt[i].p[1] << "\t" << prt[i].p[2] << endl;
            out_file->close();
            delete out_file;
        }
    }
}

void Solver::SaveOutput(int k, string file)
{
    if (k%ppw == 0)
    {
        for (int sp = 0; sp < NUM_SP; sp++)
            plasma_energy[sp] = pfc[sp].KineticEnergy(ax, ay);

        laser_energy = fdtd->Energy();

        electrostatic_energy = 0.;
        for (int i = 0; i <= mesh->MAX_Z; i++)
            electrostatic_energy += ez[i]*ez[i];
        electrostatic_energy *= 0.5*mesh->dz;

        FILE *output;
        output = fs->open_file("a+", main_dir, file.c_str());
        fprintf(output, "%06d|%16.8e|", k, ez[mesh->MAX_Z-1]);
        for (int sp = 0; sp < NUM_SP; sp++)
            fprintf(output, " %16.8e|", plasma_energy[sp]);
        fprintf(output, "%16.8e|%16.8e|%16.8e|%16.8e|\n", electrostatic_energy, laser_energy, em_flux, mm->sum(plasma_energy, NUM_SP)+electrostatic_energy+laser_energy-em_flux);
        fs->close_file(output);

        if (k==0) plasma_energy_0 = mm->sum(plasma_energy, NUM_SP);
    }
}

void Solver::SaveResults()
{
    FILE *en_out;
    en_out = fs->open_file("a+", "", "energy.txt");
    double plasma_energy_1 = mm->sum(plasma_energy, NUM_SP)+electrostatic_energy+laser_energy;
    fprintf(en_out, "%.8g\t%.8g\t%.8g\t%.8g\n", N_0, plasma_energy_0, plasma_energy_1, (plasma_energy_1-plasma_energy_0)/(0.5));
    fs->close_file(en_out);

    en_out = fs->open_file("a+", "", "energy2.txt");
    fprintf(en_out, "%.8g\t", (plasma_energy_1-plasma_energy_0)/(0.5));
    fs->close_file(en_out);
}

bool Solver::SetNotNegative(pyinput *in, string name, int *var)
{
    *var = in->GetInt(name);
    if ( *var >= 0 ) return 0;
    else {cout << name + " mustnot be negative" << endl; return 1;}
}

bool Solver::SetPositive(pyinput *in, string name, int *var)
{
    *var = in->GetInt(name);
    if ( *var > 0 ) return 0;
    else {cout << name + " must be positive" << endl; return 1;}
}

bool Solver::SetPositive(pyinput *in, string name, double *var)
{
    *var = in->GetDouble(name);
    if ( *var > 0. ) return 0;
    else {cout << name + " must be positive" << endl; return 1;}
}

bool Solver::SetFormat(pyinput *in, string name, string *var)
{
    *var = in->GetString(name);
    if (save_format == "txt" || save_format == "bin" || save_format ==  "gzip") return 0;
    else {cout << name + " must be either 'txt' or 'bin' or 'gzip'" << endl; return 1;}
}

#endif // SOLVER_H
