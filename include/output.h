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

#ifndef OUTPUT_H
#define OUTPUT_H

#include <fstream>
#include <string>

#include "pyinput.h"
#include "solver.h"

class Output {
 public:
  std::string output_directory_name;
  bool save_fields, save_concs, save_dstr;
  int save_dt, save_fields_dt, save_concs_dt, save_dstr_dt;
  std::string save_format, save_fields_format, save_concs_format,
      save_dstr_format;
  std::ofstream energy_file;

 private:
  Solver* solver;

 public:
  Output(pyinput*, std::string, Solver*, int*);
  virtual ~Output();

  // saving data
  void SaveFields(int k);    // save fields data in files
  void SaveConcs(int k);     // save concentarations data in files
  void SaveDstrFunc(int k);  // save distribution functions data in files
  void SavePrtData(int k);   // save particles data in files

  void WriteEnergy(int);

 private:
  int Init(pyinput*, std::string*);
  void CreateDirs();
  void SaveInput(std::string);
  void InitEnergyFile(std::string);
  bool SetFormat(pyinput* in, std::string name, std::string* var);
};

#endif  // OUTPUT_H
