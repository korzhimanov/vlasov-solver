/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file main.cpp
 * \brief The main source file
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 *
 * The file contains argument parser, timer, main cycle and basic logging.
 */

#ifdef _OPENMP
#include <omp.h>
#endif

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <ctime>
#include "pyinput.h"
#include "solver.h"

double invcps = 1./(double)CLOCKS_PER_SEC;
double get_time()
{
#ifdef _OPENMP
    return omp_get_wtime();
#else
    return clock()*invcps;
#endif
}

void WrongArguments()
{
    std::cout << "Wrong arguments. Type -h for help." << std::endl;
    exit(-1);
}

int main(int argc, char **argv)
{
    char input_file_name[512];
    sprintf(input_file_name, "input.py");

    char output_folder_name[512];
    time_t current_time = time(NULL);
    struct tm * utc_time = gmtime(&current_time);
    snprintf(output_folder_name, 512, "%d-%02d-%dUTC%02d:%02d:%02d", utc_time->tm_year+1900, utc_time->tm_mon+1, utc_time->tm_mday, utc_time->tm_hour, utc_time->tm_min, utc_time->tm_sec);

    if (argc > 1)
    {
        int i = 1;
        do
        {
            if (!strcmp(argv[i], "-i"))
            {
                if (argc > i+1)
                    strcpy(input_file_name, argv[++i]);
                else
                    WrongArguments();
                std::cout << "Input file name has been setted to " << argv[i] << std::endl;
                i++; continue;
            }
            if (!strcmp(argv[i], "-o"))
            {
                if (argc > i+1)
                    strcpy(output_folder_name, argv[++i]);
                else
                    WrongArguments();
                std::cout << "Output folder name has been setted to " << argv[i] << std::endl;
                i++; continue;
            }
            if (!strcmp(argv[i], "-h"))
            {
                std::cout << "Usage: vlasov [OPTION...]" << std::endl;
                std::cout << std::setw(20) << std::left << "  -i <name>" << "defines an arbitrary input file to be used" << std::endl;
                std::cout << std::setw(20) << std::left << "  -o <name>" << "defines an arbitrary output folder to be used" << std::endl;
                std::cout << std::setw(20) << std::left << "  -h" << "shows this help" << std::endl;
                exit(0);
            }
            WrongArguments();
        }
        while (i < argc);
    }

    int err = 0;

    pyinput in;
    std::cout << "Reading input file " << input_file_name << " ... " << std::endl;
    in.ReadFile(input_file_name);
    std::cout << "done!" << std::endl;

    Mesh mesh(&in, &err);
    if (err != 0) {std::cout << "Initialization failed! Returning code is " << err << std::endl; return err;}

    InitParams init_params(&in, output_folder_name, &mesh, &err);
    if (err != 0) {std::cout << "Initialization failed! Returning code is " << err << std::endl; return err;}

    Solver S(&init_params, &in, &mesh, &err);

    S.plasmas->InitDistribution();
    S.particles->InitParticles(&mesh);

    S.CreateDirs();

    S.InitOutput("output");
    S.SaveInput("input_parameters");

    double t0 = get_time();
    std::cout << std::setw(10) << "LongField" << std::setw(11) << "TransField" << std::setw(10) << "FieldGen" << std::setw(10) << "DstrFunc" << std::setw(10) << "Particles" << std::setw(10) << "SvFields" << std::setw(10) << "SvConcs" << std::setw(10) << "SvDstr" << std::setw(10) << "SvOutput" << std::setw(10) << "Step" << std::setw(20) << "Passed/ Estimated" << std::endl;

    double t1, t2, t;
    std::cout << std::fixed << std::setprecision(4);
    for (int k = 0; k < mesh.MAX_T; k++)
    {
        // LongField
        t = t1 = get_time();
        S.plasmas->CalcLongFields();
        t2 = get_time();
        std::cout << std::setw(10) << (t2-t1);

        // TransField
        t1 = t2;
        S.fdtd->CalcSources( (double)k*mesh.dt );
        S.CalcTransFields();
        t2 = get_time();
        std::cout << std::setw(11) << (t2-t1);

        // FieldGen
        t1 = t2;
        S.FieldGeneration();
        t2 = get_time();
        std::cout << std::setw(10) << (t2-t1);

        // DstrFunc
        t1 = t2;
        S.plasmas->CalcDstrFunc();
        t2 = get_time();
        std::cout << std::setw(10) << (t2-t1);

        // Particles
        t1 = t2;
        if (S.particles->particles_number > 0) S.MoveParticles();
        t2 = get_time();
        std::cout << std::setw(10) << (t2-t1);

        // SvFields
        t1 = t2;
        S.SaveFields( k );
        t2 = get_time();
        std::cout << std::setw(10) << (t2-t1);

        // SvConcs
        t1 = t2;
        S.SaveConcs( k );
        t2 = get_time();
        std::cout << std::setw(10) << (t2-t1);

        // SvDstr
        t1 = t2;
        S.SaveDstrFunc( k );
        t2 = get_time();
        std::cout << std::setw(10) << (t2-t1);

        // Output
        t1 = t2;
        //S.SaveOutput( k, "output");
        t2 = get_time();
        std::cout << std::setw(10) << (t2-t1);
        std::cout << std::setw(10) << (t2-t );
        std::cout << std::setprecision(2);
        std::cout << std::setw(9) << (t2-t0) << "/ " << (t2-t0)/(k+1)*(mesh.MAX_T) << std::endl;
        std::cout << std::setprecision(4);
    }

    //S.SaveResults();

    std::cout << "Done!" << std::endl;
    t = get_time();
    std::cout << "Full time: " << (t-t0) << " s (" << (t-t0)/mesh.MAX_T << " s per iteration)" << std::endl;

    return 0;
}
