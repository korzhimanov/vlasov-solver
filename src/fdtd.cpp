/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file fdtd.cpp
 * \brief The source file which implements the methods of the FDTD class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#include "../include/fdtd.h"

#include "include/mymath.h"

FDTD::FDTD(pyinput *in, Mesh *m, int *err)
    : mesh(m), PML(512), MAX_SIGMA(100.), SOURCE(1) {
  exl = eyl = exr = eyr = hxl = hyl = hxr = hyr = 0.;
  *err = Init(in);
  AllocMemory();
}

FDTD::~FDTD() {
  delete[] ex;
  delete[] ey;
  delete[] hx;
  delete[] hy;
  delete pulse_x;
  delete pulse_y;
  delete[] r;
  delete[] r1;
}

int FDTD::Init(pyinput *in) {
  pulse_x = new pFunc(in->GetFunc("PULSE_X"));
  pulse_y = new pFunc(in->GetFunc("PULSE_Y"));
  if (!in->SetPositive("source", &SOURCE)) return 300;

  if (!in->SetPositive("PML", &PML)) return 310;
  if (!in->SetPositive("PML_MAX_SIGMA", &MAX_SIGMA)) return 320;

  return 0;
}

void FDTD::AllocMemory() {
  ex = new double[mesh->MAX_Z + 1];
  ey = new double[mesh->MAX_Z + 1];
  hx = new double[mesh->MAX_Z];
  hy = new double[mesh->MAX_Z];
  mymath::zeros(ex, mesh->MAX_Z + 1);
  mymath::zeros(ey, mesh->MAX_Z + 1);
  mymath::zeros(hx, mesh->MAX_Z);
  mymath::zeros(hy, mesh->MAX_Z);

  r = new double[2 * PML];
  r1 = new double[2 * PML];
  CalcPMLCoeff();
}

void FDTD::Maxwell() {
  // common Yee's algorithm
  int j, p;

  // advance magnetic fields
  for (j = PML; j <= mesh->MAX_Z - PML; j++) {
    hx[j] += mesh->dt_dz * (ey[j + 1] - ey[j]);
    hy[j] -= mesh->dt_dz * (ex[j + 1] - ex[j]);
  }

  // PML with exponential time-stepping
  for (j = 0, p = 2 * PML - 1; j <= PML - 1; j++, p -= 2) {
    hx[j] = r[p] * hx[j] - r1[p] * (ey[j + 1] - ey[j]);
    hy[j] = r[p] * hy[j] + r1[p] * (ex[j + 1] - ex[j]);
  }

  for (j = mesh->MAX_Z - PML + 1, p = 1; j < mesh->MAX_Z; j++, p += 2) {
    hx[j] = r[p] * hx[j] - r1[p] * (ey[j + 1] - ey[j]);
    hy[j] = r[p] * hy[j] + r1[p] * (ex[j + 1] - ex[j]);
  }

  // incident wave correction
  hx[PML + SOURCE - 1] -= mesh->dt_dz * eyl;
  hy[PML + SOURCE - 1] += mesh->dt_dz * exl;
  hx[mesh->MAX_Z - PML - 1] += mesh->dt_dz * eyr;
  hy[mesh->MAX_Z - PML - 1] -= mesh->dt_dz * exr;

  // advance electric fields
  for (j = PML + 1; j <= mesh->MAX_Z - PML; j++) {
    ex[j] -= mesh->dt_dz * (hy[j] - hy[j - 1]);
    ey[j] += mesh->dt_dz * (hx[j] - hx[j - 1]);
  }

  // PML with exponential time-stepping
  for (j = 1, p = 2 * PML - 2; j <= PML; j++, p -= 2) {
    ex[j] = r[p] * ex[j] + r1[p] * (hy[j] - hy[j - 1]);
    ey[j] = r[p] * ey[j] - r1[p] * (hx[j] - hx[j - 1]);
  }

  for (j = mesh->MAX_Z - PML + 1, p = 0; j < mesh->MAX_Z; j++, p += 2) {
    ex[j] = r[p] * ex[j] + r1[p] * (hy[j] - hy[j - 1]);
    ey[j] = r[p] * ey[j] - r1[p] * (hx[j] - hx[j - 1]);
  }

  // incident wave correction
  ex[PML + SOURCE] += mesh->dt_dz * hyl;
  ey[PML + SOURCE] -= mesh->dt_dz * hxl;
  ex[mesh->MAX_Z - PML] += mesh->dt_dz * hyr;
  ey[mesh->MAX_Z - PML] -= mesh->dt_dz * hxr;
}

void FDTD::CalcPMLCoeff() {
  for (int j = 0; j < 2 * PML; j++) {
    double sigma = MAX_SIGMA * mymath::sqr(static_cast<double>(j + 1) /
                                           static_cast<double>(2 * PML));
    r[j] = exp(-sigma * mesh->dt);
    r1[j] = (r[j] - 1.) / (sigma * mesh->dz);
  }
}
