/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file particle.cpp
 * \brief The source file which implements the methods of the Particle class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#include "include/mesh.h"

Mesh::Mesh()
    : ppw(16),
      MAX_Z(32),
      MAX_T(32),
      dz(2. * M_PI / 16),
      dt(2. * M_PI / 16),
      dt_dz(1.) {}

Mesh::Mesh(pyinput *in, int *err)
    : ppw(16),
      MAX_Z(32),
      MAX_T(32),
      dz(2. * M_PI / 16),
      dt(2. * M_PI / 16),
      dt_dz(1.) {
  *err = Init(in);
}

int Mesh::Init(pyinput *in) {
  if (!in->SetPositive("ppw", &ppw)) return 100;
  if (!in->SetPositive("MAX_Z", &MAX_Z)) return 110;
  if (!in->SetPositive("MAX_T", &MAX_T)) return 120;
  if (!in->SetPositive("dz", &dz)) return 130;
  if (!in->SetPositive("dt", &dt)) return 140;
  dt_dz = dt / dz;
  return 0;
}
