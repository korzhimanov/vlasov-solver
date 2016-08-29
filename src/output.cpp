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

#include "include/output.h"

#include <iomanip>
#include <iostream>
#include <sstream>

#include "include/errors.h"
#include "include/file_saving.h"
#include "include/mymath.h"

Output::Output(const PyInput &in, std::string dn, Solver *s, int &err)
    : save_format("gzip"),
      save_fields_format(""),
      save_concs_format(""),
      save_dstr_format(""),
      solver(s) {
  output_directory_name = dn + "/";

  in.Set("save_fields", save_fields, err);
  if (err == VAR_NOT_FOUND) {
    std::cout << "Set save_fields to 'true'." << std::endl;
    save_fields = true;
    err = 0;
  }
  in.Set("save_concs", save_concs, err);
  if (err == VAR_NOT_FOUND) {
    std::cout << "Set save_concs to 'true'." << std::endl;
    save_concs = true;
    err = 0;
  }
  in.Set("save_dstr", save_dstr, err);
  if (err == VAR_NOT_FOUND) {
    std::cout << "Set save_dstr to 'false'." << std::endl;
    save_dstr = false;
    err = 0;
  }

  if ((save_fields || save_concs || save_dstr) == true) {
    in.SetPositive("save_dt", save_dt, err);
    if (err == VAR_NOT_FOUND) {
      std::cout << "Set save_dt to '1'." << std::endl;
      save_dt = 1;
      err = 0;
    }
    SetFormat(in, "save_format", save_format, err);
    if (err == VAR_NOT_FOUND) {
      std::cout << "Set save_format to 'txt'." << std::endl;
      save_format = "txt";
      err = 0;
    }
    if (save_format == "") {
      save_format = "txt";
    }

    if (save_fields == true) {
      in.SetNotNegative("save_fields_dt", save_fields_dt, err);
      if (err == VAR_NOT_FOUND) {
        std::cout << "Set save_fields_dt to " << save_dt << "." << std::endl;
        save_fields_dt = save_dt;
        err = 0;
      }
      if (save_fields_dt == 0) save_fields_dt = save_dt;
      SetFormat(in, "save_fields_format", save_fields_format, err);
      if (err == VAR_NOT_FOUND) {
        std::cout << "Set save_fields_format to " << save_format << "."
                  << std::endl;
        save_fields_format = save_format;
        err = 0;
      }
      if (save_fields_format == "") save_fields_format = save_format;
    }
    if (save_concs == true) {
      in.SetNotNegative("save_concs_dt", save_concs_dt, err);
      if (err == VAR_NOT_FOUND) {
        std::cout << "Set save_concs_dt to " << save_dt << "." << std::endl;
        save_concs_dt = save_dt;
        err = 0;
      }
      if (save_concs_dt == 0) save_concs_dt = save_dt;
      SetFormat(in, "save_concs_format", save_concs_format, err);
      if (err == VAR_NOT_FOUND) {
        std::cout << "Set save_concs_format to " << save_format << "."
                  << std::endl;
        save_concs_format = save_format;
        err = 0;
      }
      if (save_concs_format == "") save_concs_format = save_format;
    }
    if (save_dstr == true) {
      in.SetNotNegative("save_dstr_dt", save_dstr_dt, err);
      if (err == VAR_NOT_FOUND) {
        std::cout << "Set save_dstr_dt to " << save_dt << "." << std::endl;
        save_dstr_dt = save_dt;
        err = 0;
      }
      if (save_dstr_dt == 0) save_dstr_dt = save_dt;
      SetFormat(in, "save_dstr_format", save_dstr_format, err);
      if (err == VAR_NOT_FOUND) {
        std::cout << "Set save_dstr_format to " << save_format << "."
                  << std::endl;
        save_dstr_format = save_format;
        err = 0;
      }
      if (save_dstr_format == "") save_dstr_format = save_format;
    }
  }

  CreateDirs();
  SaveInput("initial_parameters.txt");
  InitEnergyFile("energy.txt");
}

Output::~Output() {
  if (energy_file) energy_file.close();
}

void Output::SetFormat(const PyInput &in, std::string name, std::string &var,
                       int &err) {
  in.Set(name, var, err);
  if (err != 0) return;

  if (var != "" && var != "txt" && var != "bin" && var != "gzip") {
    std::cout << name + " must be either '', 'txt', 'bin' or 'gzip'"
              << std::endl;
    err = WRONG_OUTPUT_FORMAT;
  }
}

void Output::CreateDirs() {
  // creating main directory for output
  { filesaving::create_dir(output_directory_name); }

  // creating directory for test particles
  if (solver->particles->particles_number > 0) {
    filesaving::create_dir(output_directory_name + "particles/");
  }

  // creating directories for storing distribution functions in files
  if (save_dstr == true) {
    for (int sp = 0; sp < solver->plasmas->species_number; sp++) {
      std::ostringstream ss;
      ss << output_directory_name << "dstr_func" << sp << "/";
      filesaving::create_dir(ss.str());
    }
  }
}

void Output::SaveInput(std::string fn) {
  // storing information about input parameters
  std::ofstream fs((output_directory_name + fn).c_str(), std::ios_base::app);
  fs << "Relativistic Vlasov-Maxwell solver. Developed by Artem Korzhimanov. "
        "mailto: korzhimanov.artem@gmail.com\n\n";

  fs << "GENERAL PARAMETERS\n\n";
  fs << "\tSpace cells per wavelength = " << solver->mesh->ppw << "\n";
  fs << "\tSpace step = " << .5 * M_1_PI *solver->mesh->dz << " Wavelengths\n";
  fs << "\tTotal space cells = " << solver->mesh->MAX_Z << "\n";
  fs << "\tTotal space size = " << static_cast<double>(solver->mesh->MAX_Z) /
                                       solver->mesh->ppw << "Wavelengths\n";
  fs << "\tTime step = " << .5 * M_1_PI *solver->mesh->dt << " Waveperiods\n";
  fs << "\tTotal time steps = " << solver->mesh->MAX_T << "\n";
  fs << "\tTotal running time = "
     << .5 * M_1_PI *solver->mesh->MAX_T *solver->mesh->dt << " Waveperiods\n";

  fs << "\nPLASMA PARAMETERS\n\n";
  for (int sp = 0; sp < solver->plasmas->species_number; sp++) {
    fs << "\n\tSpecies " << sp << "\n\n";
    solver->plasmas->pfc[sp].SaveInput(fs);
  }

  fs << "\nTEST PARTICLES PARAMETERS\n\n";
  fs << "\tNumber of particles = " << solver->particles->particles_number
     << "\n";
  if (solver->particles->particles_number > 0) {
    fs << "\tMass of particles = " << solver->particles->mass << "\n";
    fs << "\tCharge of particles = " << solver->particles->charge << "\n";
    fs << "\tParticles are equidistantly deposited from "
       << static_cast<double>(solver->particles->start_point) /
              solver->mesh->ppw;
    fs << " to "
       << static_cast<double>(solver->particles->start_point +
                              solver->particles->interval *
                                  solver->particles->particles_number) /
              solver->mesh->ppw << " Wavelengths\n";
  }

  fs << "\nOUTPUT PARAMETERS\n\n";
  fs << "\tFields are ";
  if (save_fields == true) {
    fs << "being saved every " << save_fields_dt << " steps";
    fs << " or every " << static_cast<double>(save_fields_dt) /
                              solver->mesh->ppw << " Waveperiods in files";
    fs << " in " << save_fields_format << "-format\n";
  } else
    fs << "not being saved\n";
  fs << "\tConcentration distributions are ";
  if (save_concs == true) {
    fs << "being saved every " << save_concs_dt << " steps";
    fs << " or every " << static_cast<double>(save_concs_dt) / solver->mesh->ppw
       << " Waveperiods in files";
    fs << " in " << save_concs_format << "-format\n";
  } else
    fs << "not being saved\n";
  fs << "\tDistribution functions are ";
  if (save_dstr == true) {
    fs << "being saved every " << save_dstr_dt << " steps";
    fs << " or every " << static_cast<double>(save_dstr_dt) / solver->mesh->ppw
       << " Waveperiods in files";
    fs << " in " << save_dstr_format << "-format";
  } else
    fs << "not being saved";

  fs << std::endl;
}

void Output::SaveFields(int k) {
  if (save_fields == true && k % save_fields_dt == 0) {
    if (save_fields_format == "txt") {
      filesaving::save_file_1D_txt(solver->fdtd->ex, solver->mesh->MAX_Z,
                                   output_directory_name + "ex.txt");
      filesaving::save_file_1D_txt(solver->fdtd->ey, solver->mesh->MAX_Z,
                                   output_directory_name + "ey.txt");
      filesaving::save_file_1D_txt(solver->fdtd->hx, solver->mesh->MAX_Z,
                                   output_directory_name + "hx.txt");
      filesaving::save_file_1D_txt(solver->fdtd->hy, solver->mesh->MAX_Z,
                                   output_directory_name + "hy.txt");
      filesaving::save_file_1D_txt(solver->plasmas->ez, solver->mesh->MAX_Z,
                                   output_directory_name + "ez.txt");
      filesaving::save_file_1D_txt(solver->plasmas->a2, solver->mesh->MAX_Z,
                                   output_directory_name + "vecpot2.txt");
    }

    if (save_fields_format == "bin") {
      filesaving::save_file_1D_bin(solver->fdtd->ex, solver->mesh->MAX_Z,
                                   output_directory_name + "ex.bin");
      filesaving::save_file_1D_bin(solver->fdtd->ey, solver->mesh->MAX_Z,
                                   output_directory_name + "ey.bin");
      filesaving::save_file_1D_bin(solver->fdtd->hx, solver->mesh->MAX_Z,
                                   output_directory_name + "hx.bin");
      filesaving::save_file_1D_bin(solver->fdtd->hy, solver->mesh->MAX_Z,
                                   output_directory_name + "hy.bin");
      filesaving::save_file_1D_bin(solver->plasmas->ez, solver->mesh->MAX_Z,
                                   output_directory_name + "ez.bin");
      filesaving::save_file_1D_bin(solver->plasmas->a2, solver->mesh->MAX_Z,
                                   output_directory_name + "vecpot2.bin");
    }

    if (save_fields_format == "gzip") {
      filesaving::save_file_1D_gzip(solver->fdtd->ex, solver->mesh->MAX_Z,
                                    output_directory_name + "ex.gz");
      filesaving::save_file_1D_gzip(solver->fdtd->ey, solver->mesh->MAX_Z,
                                    output_directory_name + "ey.gz");
      filesaving::save_file_1D_gzip(solver->fdtd->hx, solver->mesh->MAX_Z,
                                    output_directory_name + "hx.gz");
      filesaving::save_file_1D_gzip(solver->fdtd->hy, solver->mesh->MAX_Z,
                                    output_directory_name + "hy.gz");
      filesaving::save_file_1D_gzip(solver->plasmas->ez, solver->mesh->MAX_Z,
                                    output_directory_name + "ez.gz");
      filesaving::save_file_1D_gzip(solver->plasmas->a2, solver->mesh->MAX_Z,
                                    output_directory_name + "vecpot2.gz");
    }
  }
}

void Output::SaveConcs(int k) {
  if (save_concs == true && k % save_concs_dt == 0) {
    for (int sp = 0; sp < solver->plasmas->species_number; sp++) {
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

void Output::SaveDstrFunc(int k) {
  if (save_dstr == true && k % save_dstr_dt == 0) {
    for (int sp = 0; sp < solver->plasmas->species_number; sp++) {
      std::stringstream ss;
      ss << output_directory_name << "dstr_func" << sp << "/data"
         << std::setfill('0') << std::setw(8) << k;
      if (save_dstr_format == "txt")
        solver->plasmas->pfc[sp].SaveDstrFunctionTxt(ss.str());

      if (save_dstr_format == "bin")
        solver->plasmas->pfc[sp].SaveDstrFunctionBin(ss.str());

      if (save_dstr_format == "gzip")
        solver->plasmas->pfc[sp].SaveDstrFunctionGZip(ss.str());
    }
  }
}

void Output::SavePrtData(int k) {
  if (solver->particles->particles_number > 0) {
    for (int i = 0; i < solver->particles->particles_number; i++) {
      std::ofstream *out_file;
      std::stringstream file_path;
      file_path << output_directory_name << "particles/particle"
                << std::setfill('0') << std::setw(6) << i << ".txt";
      out_file = new std::ofstream(file_path.str().c_str(), std::ios_base::app);
      *out_file << solver->particles->prt[i].r[0] << "\t"
                << solver->particles->prt[i].r[1] << "\t"
                << solver->particles->prt[i].r[2] << "\t"
                << solver->particles->prt[i].p[0] << "\t"
                << solver->particles->prt[i].p[1] << "\t"
                << solver->particles->prt[i].p[2] << std::endl;
      out_file->close();
      delete out_file;
    }
  }
}

void Output::InitEnergyFile(std::string fn) {
  energy_file.open((output_directory_name + fn).c_str(), std::ios_base::app);
  energy_file << std::right << std::scientific << std::setprecision(4);
  energy_file << std::setw(8) << "Step";
  energy_file << std::setw(16) << "Neutrality";
  for (int sp = 0; sp < solver->plasmas->species_number; sp++)
    energy_file << std::setw(9) << sp << " specie";
  energy_file << std::setw(16) << "Electrostatic" << std::setw(16) << "Laser"
              << std::setw(16) << "EM Flux" << std::setw(16) << "Total"
              << std::setw(16) << "Balance" << std::endl;
}

void Output::WriteEnergy(int k) {
  // initialize variables needed for function on the very first call
  static double em_flux = 0., electrostatic_energy = 0., laser_energy = 0.,
                *plasma_energy, plasma_energy_0 = 0., total_energy = 0.;

  plasma_energy = new double[solver->plasmas->species_number];
  em_flux += solver->fdtd->FluxIn();

  if (k % (solver->mesh->ppw / 16) == 0) {
    for (int sp = 0; sp < solver->plasmas->species_number; sp++)
      plasma_energy[sp] = solver->plasmas->pfc[sp].KineticEnergy(
          solver->plasmas->ax, solver->plasmas->ay);

    if (k == 0)
      plasma_energy_0 =
          mymath::sum(plasma_energy, solver->plasmas->species_number);

    laser_energy = solver->fdtd->Energy();

    electrostatic_energy = 0.;
    for (int i = 0; i <= solver->mesh->MAX_Z; i++)
      electrostatic_energy += solver->plasmas->ez[i] * solver->plasmas->ez[i];
    electrostatic_energy *= 0.5 * solver->mesh->dz;

    energy_file << std::setw(8) << k;
    energy_file << std::setw(16)
                << solver->plasmas->ez[solver->mesh->MAX_Z - 1];
    for (int sp = 0; sp < solver->plasmas->species_number; sp++)
      energy_file << std::setw(16) << plasma_energy[sp];

    total_energy = mymath::sum(plasma_energy, solver->plasmas->species_number) +
                   electrostatic_energy + laser_energy;
    energy_file << std::setw(16) << electrostatic_energy << std::setw(16)
                << laser_energy << std::setw(16) << em_flux << std::setw(16)
                << total_energy << std::setw(16)
                << total_energy - em_flux - plasma_energy_0 << std::endl;
  }
  delete[] plasma_energy;
}
