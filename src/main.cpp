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

#include <cstring>
#include <ctime>

#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "include/errors.h"
#include "include/mesh.h"
#include "include/output.h"
#include "include/pyinput.h"
#include "include/solver.h"

const double invcps = 1. / static_cast<double>(CLOCKS_PER_SEC);
double get_time() {
#ifdef _OPENMP
  return omp_get_wtime();
#else
  return clock() * invcps;
#endif
}

void WrongArguments() {
  std::cerr << "Wrong arguments. Type -h for help." << std::endl;
  exit(-1);
}

int main(int argc, char **argv) {
  std::stringstream init_file_name("init/default.py");

  std::stringstream output_folder_name;
  time_t current_time = time(NULL);
  struct tm *utc_time = gmtime(&current_time);

  output_folder_name << utc_time->tm_year + 1900 << "-";
  output_folder_name << std::setfill('0');
  output_folder_name << std::setw(2) << utc_time->tm_mon + 1 << "-"
                     << std::setw(2) << utc_time->tm_mday << "UTC"
                     << std::setw(2) << utc_time->tm_hour << ":" << std::setw(2)
                     << utc_time->tm_min << ":" << std::setw(2)
                     << utc_time->tm_sec;

  if (argc > 1) {
    int i = 1;
    do {
      if (!strcmp(argv[i], "-i")) {
        if (argc > i + 1) {
          init_file_name.str(std::string());
          init_file_name.clear();
          init_file_name << argv[++i];
        } else
          WrongArguments();
        std::cout << "Init-file name has been setted to " << argv[i]
                  << std::endl;
        i++;
        continue;
      }
      if (!strcmp(argv[i], "-o")) {
        if (argc > i + 1) {
          output_folder_name.str(std::string());
          output_folder_name.clear();
          output_folder_name << argv[++i];
        } else
          WrongArguments();
        std::cout << "Output folder name has been setted to " << argv[i]
                  << std::endl;
        i++;
        continue;
      }
      if (!strcmp(argv[i], "-h")) {
        std::cout << "Usage: vlasov [OPTION...]" << std::endl;
        std::cout << std::setw(20) << std::left << "  -i <name>"
                  << "defines an arbitrary input file to be used" << std::endl;
        std::cout << std::setw(20) << std::left << "  -o <name>"
                  << "defines an arbitrary output folder to be used"
                  << std::endl;
        std::cout << std::setw(20) << std::left << "  -h"
                  << "shows this help" << std::endl;
        exit(0);
      }
      WrongArguments();
    } while (i < argc);
  }

  int err = 0;

  PyInput in(err);
  if (err) exit(err);
  std::cout << "Reading init-file " << init_file_name.str() << " ..."
            << std::endl;
  in.ReadFile(init_file_name.str(), err);
  if (err) {
    std::cerr << "failed!\n Terminated!" << std::endl;
    exit(err);
  }
  std::cout << "done!" << std::endl;

  std::cout << "Initalizing mesh ..." << std::endl;
  Mesh mesh(in, err);
  if (err) {
    std::cerr << "failed!\n Terminated!" << std::endl;
    exit(err);
  }
  std::cout << "done!" << std::endl;

  std::cout << "Initializing solver ..." << std::endl;
  Solver S(in, &mesh, err);
  if (err) {
    std::cerr << "failed!\n Terminated!" << std::endl;
    exit(err);
  }
  std::cout << "done!" << std::endl;

  std::cout << "Initializing output ..." << std::endl;
  Output output(in, output_folder_name.str(), &S, err);
  if (err) {
    std::cerr << "failed!\n Terminated!" << std::endl;
    exit(err);
  }
  std::cout << "Initialization of output is done!" << std::endl;

  S.plasmas->InitDistribution();
  S.particles->InitParticles(&mesh);

  double t0 = get_time();
  std::cout << std::setw(10) << "LongField" << std::setw(11) << "TransField"
            << std::setw(10) << "FieldGen" << std::setw(10) << "DstrFunc"
            << std::setw(10) << "Particles" << std::setw(10) << "SvFields"
            << std::setw(10) << "SvConcs" << std::setw(10) << "SvDstr"
            << std::setw(10) << "SvPrtData " << std::setw(10) << " Energy "
            << std::setw(10) << "Step" << std::setw(20) << "Passed/ Estimated"
            << std::endl;

  double t;
  std::cout << std::fixed << std::setprecision(4);
  for (int k = 0; k < mesh.MAX_T; k++) {
    double t1, t2;

    // LongField
    t = t1 = get_time();
    S.plasmas->CalcLongFields();
    t2 = get_time();
    std::cout << std::setw(10) << (t2 - t1);

    // TransField
    t1 = t2;
    S.fdtd->CalcSources(static_cast<double>(k) * mesh.dt);
    S.CalcTransFields();
    t2 = get_time();
    std::cout << std::setw(11) << (t2 - t1);

    // FieldGen
    t1 = t2;
    S.FieldGeneration();
    t2 = get_time();
    std::cout << std::setw(10) << (t2 - t1);

    // DstrFunc
    t1 = t2;
    S.plasmas->CalcDstrFunc();
    t2 = get_time();
    std::cout << std::setw(10) << (t2 - t1);

    // Particles
    t1 = t2;
    if (S.particles->particles_number > 0) S.MoveParticles();
    t2 = get_time();
    std::cout << std::setw(10) << (t2 - t1);

    // SvFields
    t1 = t2;
    output.SaveFields(k);
    t2 = get_time();
    std::cout << std::setw(10) << (t2 - t1);

    // SvConcs
    t1 = t2;
    output.SaveConcs(k);
    t2 = get_time();
    std::cout << std::setw(10) << (t2 - t1);

    // SvDstr
    t1 = t2;
    output.SaveDstrFunc(k);
    t2 = get_time();
    std::cout << std::setw(10) << (t2 - t1);

    // SvPrtData
    t1 = t2;
    output.SavePrtData(k);
    t2 = get_time();
    std::cout << std::setw(10) << (t2 - t1);

    // Output
    t1 = t2;
    output.WriteEnergy(k);
    t2 = get_time();
    std::cout << std::setw(10) << (t2 - t1);
    std::cout << std::setw(10) << (t2 - t);
    std::cout << std::setprecision(2);
    std::cout << std::setw(9) << (t2 - t0) << "/ "
              << (t2 - t0) / (k + 1) * (mesh.MAX_T) << std::endl;
    std::cout << std::setprecision(4);
  }

  std::cout << "Done!" << std::endl;
  t = get_time();
  std::cout << "Full time: " << (t - t0) << " s (" << (t - t0) / mesh.MAX_T
            << " s per iteration)" << std::endl;

  return 0;
}
