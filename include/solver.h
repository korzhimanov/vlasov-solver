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
 * \bug The methods of the Solver class are implemented in the header due to
 * some problems with compilation.
 */

#ifndef INCLUDE_SOLVER_H_
#define INCLUDE_SOLVER_H_

#include "include/fdtd.h"
#include "include/mesh.h"
#include "include/plasmas.h"
#include "include/pyinput.h"
#include "include/testparticles.h"

/**
 * \class Solver
 * The main class where actually the main sub-steps of each iteration are
 * defined
 */
class Solver {
 public:
  Mesh *mesh;
  // variables containing data
  Plasmas *plasmas;
  FDTD *fdtd;
  TestParticles *particles;

 private:
  // auxiliary variables
  double halfn0dz;

 public:
  // ------constructors---------------------------------------------------
  Solver(pyinput *, Mesh *, int *);
  // ------destructor-----------------------------------------------------
  virtual ~Solver();

  // ------initialisation-------------------------------------------------
  int Init(pyinput *, int *);

  // ------calculations---------------------------------------------------
  inline void CalcTransFields() {
    // evaluates vector potential
    for (int i = 0; i <= mesh->MAX_Z; i++) {
      plasmas->ax[i] -= fdtd->ex[i] * mesh->dt;
      plasmas->ay[i] -= fdtd->ey[i] * mesh->dt;
      plasmas->a2[i] =
          plasmas->ax[i] * plasmas->ax[i] + plasmas->ay[i] * plasmas->ay[i];
    }

    // evaluating electric and magnetic fields
    fdtd->Maxwell();
  }

  void MoveParticles();

  inline void FieldGeneration() {
    // field generation by plasma currents
    for (int sp = 0; sp < plasmas->species_number; sp++)
      plasmas->pfc[sp].CalcCurrent(fdtd, plasmas->ax, plasmas->ay, plasmas->a2);

    // field generation by fixed ions (only in a boosted frame)
    if (halfn0dz > 0.) {
      fdtd->ey[0] -= halfn0dz * plasmas->fixed_ions_conc[0];
      for (int i = 1; i < mesh->MAX_Z; i++)
        fdtd->ey[i] -= halfn0dz * (plasmas->fixed_ions_conc[i] +
                                   plasmas->fixed_ions_conc[i - 1]);
      fdtd->ey[mesh->MAX_Z] -=
          halfn0dz * plasmas->fixed_ions_conc[mesh->MAX_Z - 1];
    }
  }

 private:
};

#endif  // INCLUDE_SOLVER_H_
