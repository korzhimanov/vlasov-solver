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

#include "solver.h"
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <ctime>

using namespace std;

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
    cout << "Wrong arguments. Type -h for help." << endl; exit(-1);
}

int main(int argc, char **argv)
{
    char input_file_name[512];
    sprintf(input_file_name, "input.py");

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
                cout << "Input file name has been setted to " << argv[i] << endl;
                i++; continue;
            }
            if (!strcmp(argv[i], "-h"))
            {
                cout << "Usage: vlasov [OPTION...]" << endl;
                cout << setw(20) << left << "  -i <name>" << "defines an arbitrary input file to be used" << endl;
                cout << setw(20) << left << "  -h" << "shows this help" << endl;
                exit(0);
            }
            WrongArguments();
        }
        while (i < argc);
    }

    Solver S(input_file_name);

    S.InitFields();
    S.InitPlasma();
    S.InitTestParticles();

    S.CreateDirs();

    S.SaveInput("input_parameters");
    S.InitOutput("output");

    double t0 = get_time();
    cout << setw(10) << "LongField" << setw(11) << "TransField" << setw(10) << "FieldGen" << setw(10) << "DstrFunc" << setw(10) << "Particles" << setw(10) << "SvFields" << setw(10) << "SvConcs" << setw(10) << "SvDstr" << setw(10) << "SvOutput" << setw(10) << "Step" << setw(20) << "Passed/ Estimated" << endl;

    double t1, t2, t;
    cout << fixed << setprecision(4);
    for (int k = 0; k < S.mesh->MAX_T; k++)
    {
        // LongField
        t = t1 = get_time();
        S.CalcLongFields();
        t2 = get_time();
        cout << setw(10) << (t2-t1);

        // TransField
        t1 = t2;
        S.CalcSources( (double)k*S.mesh->dt );
        S.CalcTransFields();
        t2 = get_time();
        cout << setw(11) << (t2-t1);

        // FieldGen
        t1 = t2;
        S.FieldGeneration();
        t2 = get_time();
        cout << setw(10) << (t2-t1);

        // DstrFunc
        t1 = t2;
        S.CalcDstrFunc();
        t2 = get_time();
        cout << setw(10) << (t2-t1);

        // Particles
        t1 = t2;
        if (S.NUM_PRT > 0) S.MoveParticles();
        t2 = get_time();
        cout << setw(10) << (t2-t1);

        // SvFields
        t1 = t2;
        S.SaveFields( k );
        t2 = get_time();
        cout << setw(10) << (t2-t1);

        // SvConcs
        t1 = t2;
        S.SaveConcs( k );
        t2 = get_time();
        cout << setw(10) << (t2-t1);

        // SvDstr
        t1 = t2;
        S.SaveDstrFunc( k );
        t2 = get_time();
        cout << setw(10) << (t2-t1);

        // Output
        t1 = t2;
        //S.SaveOutput( k, "output");
        t2 = get_time();
        cout << setw(10) << (t2-t1);
        cout << setw(10) << (t2-t );
        cout << setprecision(2);
        cout << setw(9) << (t2-t0) << "/ " << (t2-t0)/(k+1)*(S.mesh->MAX_T) << endl;
        cout << setprecision(4);
    }

    //S.SaveResults();

    cout << "Done!" << endl;
    t = get_time();
    cout << "Full time: " << (t-t0) << " s (" << (t-t0)/S.mesh->MAX_T << " s per iteration)" << endl;

    return 0;
}
