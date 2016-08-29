/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file output.h
 * \brief The file defines Output class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 *
 * The file contains argument parser, timer, main cycle and basic logging.
 */

#ifndef INCLUDE_OUTPUT_H_
#define INCLUDE_OUTPUT_H_

#include <fstream>
#include <string>

#include "include/pyinput.h"
#include "include/solver.h"

class Output {
 public:
  std::string output_directory_name;
  bool save_fields, save_concs, save_dstr;
  int save_dt, save_fields_dt, save_concs_dt, save_dstr_dt;
  std::string save_format, save_fields_format, save_concs_format,
      save_dstr_format;
  std::ofstream energy_file;

 private:
  Solver *solver;

 public:
  Output(const PyInput &, std::string output_directory_name, Solver *,
         int &err);
  virtual ~Output();

  // saving data
  void SaveFields(int k);    // save fields data in files
  void SaveConcs(int k);     // save concentarations data in files
  void SaveDstrFunc(int k);  // save distribution functions data in files
  void SavePrtData(int k);   // save particles data in files

  void WriteEnergy(int);

 private:
  void CreateDirs();
  void SaveInput(std::string);
  void InitEnergyFile(std::string);
  void SetFormat(const PyInput &, std::string save_format, std::string &var,
                 int &err);
};

#endif  // INCLUDE_OUTPUT_H_
