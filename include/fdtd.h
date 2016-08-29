/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file fdtd.h
 * \brief The header file which defines FDTD class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#ifndef INCLUDE_FDTD_H_
#define INCLUDE_FDTD_H_

#include "include/mesh.h"
#include "include/pfunc.h"
#include "include/pyinput.h"

/**
 * /class FDTD
 * The class implements the FDTD solver
*/
class FDTD {
 public:
  double *ex, *ey, *hx, *hy;  // fields
 private:
  double exl, eyl, exr, eyr, hxl, hyl, hxr, hyr;  // field sources on boundaries

 private:
  Mesh *mesh;
  int PML;           // number of PML-steps
  double MAX_SIGMA;  // maximal absorption coefficient in PML
  double *r, *r1;    // absorption coefficients in PML
  int SOURCE;        // the position of a source relative to PML layer
  pFunc *pulse_x, *pulse_y;

 public:
  FDTD(pyinput *, Mesh *, int *err);
  virtual ~FDTD();

 private:
  // -----initialization-------------------------------------------------
  // initializes all but PML parameters, allocates memory for fields and puts
  // all of them equal to zero
  int Init(pyinput *);
  // allocates memory for fields and puts all of them equal to zero
  void AllocMemory();

 public:
  // -----solver---------------------------------------------------------
  // solves Maxwell equations
  void Maxwell();

  /**
   * \todo Change .5*mesh->dz and .5*mesh->dt to constants
   */
  inline void CalcSources(double t) {
    exl = pulse_x->call(t - .5 * mesh->dz);
    eyl = pulse_y->call(t - .5 * mesh->dz);
    hxl = -pulse_y->call(t + .5 * mesh->dt);
    hyl = pulse_x->call(t + .5 * mesh->dt);
  }

  // -----informative functions------------------------------------------
  // calculates flux in x-direction at 'cell'
  inline double Flux(int cell) {
    return 0.5 * (ex[cell] * (hy[cell] + hy[cell - 1]) -
                  ey[cell] * (hx[cell] + hx[cell - 1])) *
           mesh->dt;
  }
  // calculates flux incoming into box through boundaries
  inline double FluxIn() {
    return Flux(PML + SOURCE + 2) - Flux(mesh->MAX_Z - PML - 3);
  }
  // calculates flux outcoming out of box through boundaries
  inline double FluxOut() { return -FluxIn(); }
  // calculates full electromagnetic energy between 'l_cell' and 'r_cell'
  // by default calculates energy consisted in the region between PML regions
  inline double Energy(int l_cell, int r_cell) {
    double tmp = 0.;
    for (int i = l_cell; i < r_cell; i++)
      tmp += hx[i] * hx[i] + hy[i] * hy[i] + ex[i] * ex[i] + ey[i] * ey[i];
    tmp += 0.5 * (ex[r_cell] * ex[r_cell] + ey[r_cell] * ey[r_cell] -
                  ex[l_cell] * ex[l_cell] - ey[l_cell] * ey[l_cell]);
    tmp *= 0.5 * mesh->dz;
    return tmp;
  }
  inline double Energy() {
    return Energy(PML + SOURCE + 2, mesh->MAX_Z - PML - 3);
  }

 private:
  // -----miscellaneous--------------------------------------------------
  // calculates PML absorption coefficients
  void CalcPMLCoeff();
};

#endif  // INCLUDE_FDTD_H_
