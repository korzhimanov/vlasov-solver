/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file plasmas.h
 * \brief The header file which defines Plasmas class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#ifndef PLASMAS_H
#define PLASMAS_H

#include "include/mesh.h"
#include "include/mymath.h"
#include "include/pfc.h"
#include "include/pfunc.h"
#include "include/pyinput.h"

class Plasmas {
 public:
  int species_number;
  double critical_concentration,  // critical concentration in a boosted frame
                                  // (useful for simulation of oblique
                                  // irradiation)
      mean_initial_momentum;  // initial mean momentum along y-axis (useful for
                              // simulation of oblique irradiation in boosted
                              // frame)
  pFunc *fixed_ions_profile;
  double *fixed_ions_conc;  // concentration of fixed ions
  PFC *pfc;                 // PFC - class
  double *ax, *ay, *a2;     // vector potential
  double *ez;               // longitudinal field

 public:
  Plasmas(pyinput *, Mesh *, int *);
  virtual ~Plasmas();
  int Init(pyinput *);      // initializes plasma distribution
  void AllocMemory();       // allocates memory
  void InitDistribution();  // initializes distribution of plasmas
  void CalcLongFields();
  void CalcDstrFunc();

 private:
  Mesh *mesh;
};

/**
 * \todo Change 0.5*critical_concentration*mesh->dz to a single constant
 */
inline void Plasmas::CalcLongFields() {
  // calculates longitudinal field
  mymath::zeros(ez, mesh->MAX_Z + 1);
  for (int sp = 0; sp < species_number; sp++) pfc[sp].CalcLongitudinalField(ez);

  ez[0] += 0.5 * critical_concentration * mesh->dz * fixed_ions_conc[0];
  for (int i = 1; i < mesh->MAX_Z; i++)
    ez[i] += ez[i - 1] +
             0.5 * critical_concentration * mesh->dz *
                 (fixed_ions_conc[i] + fixed_ions_conc[i - 1]);
  ez[mesh->MAX_Z] += 0.5 * critical_concentration * mesh->dz *
                     fixed_ions_conc[mesh->MAX_Z - 1];
}

inline void Plasmas::CalcDstrFunc() {
  for (int sp = 0; sp < species_number; sp++) pfc[sp].MakeStep(ez, ax, ay);
}

#endif  // PLASMAS_H
