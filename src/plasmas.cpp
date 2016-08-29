/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file plasmas.cpp
 * \brief The source file which implements the methods of the Plasmas class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#include "include/plasmas.h"

#include <iostream>

#include "include/errors.h"

Plasmas::Plasmas(const PyInput &in, Mesh *m, const double &theta, int &err)
    : species_number(0),
      critical_concentration(1.),
      mean_initial_momentum(0.),
      mesh(m) {
  fixed_ions_profile = new pFunc(in.GetFunc("FIXED_IONS_PROFILE", err));

  in.SetNotNegative("NUM_SP", species_number, err);

  critical_concentration = 1. / (cos(theta) * cos(theta) * cos(theta));
  mean_initial_momentum = tan(theta);

  pfc = new PFC[species_number];
  for (int sp = 0; sp < species_number; sp++) {
    pfc[sp].Init(sp, in, mesh, &critical_concentration, &mean_initial_momentum,
                 err);
  }
  AllocMemory();
}

Plasmas::~Plasmas() {
  if (fixed_ions_profile) delete fixed_ions_profile;
  delete[] fixed_ions_conc;
  delete[] pfc;
  delete[] ax;
  delete[] ay;
  delete[] a2;
}

void Plasmas::AllocMemory() {
  fixed_ions_conc = new double[mesh->MAX_Z];

  for (int sp = 0; sp < species_number; sp++) pfc[sp].AllocMemory();

  ax = new double[mesh->MAX_Z + 1];
  ay = new double[mesh->MAX_Z + 1];
  a2 = new double[mesh->MAX_Z + 1];
  ez = new double[mesh->MAX_Z + 1];
}

void Plasmas::InitDistribution() {
  for (int i = 0; i < mesh->MAX_Z; i++)
    fixed_ions_conc[i] = fixed_ions_profile->call(i);
  delete fixed_ions_profile;  // not needed anymore

  for (int sp = 0; sp < species_number; sp++) pfc[sp].SetDistribution();

  mymath::zeros(ax, mesh->MAX_Z + 1);
  mymath::zeros(ay, mesh->MAX_Z + 1);
  mymath::zeros(a2, mesh->MAX_Z + 1);
  mymath::zeros(ez, mesh->MAX_Z + 1);
}
