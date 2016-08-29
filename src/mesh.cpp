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

Mesh::Mesh(const PyInput &in, int &err)
    : ppw(16),
      MAX_Z(32),
      MAX_T(32),
      dz(2. * M_PI / 16),
      dt(2. * M_PI / 16),
      dt_dz(1.) {
  in.SetPositive("ppw", ppw, err);
  in.SetPositive("MAX_Z", MAX_Z, err);
  in.SetPositive("MAX_T", MAX_T, err);
  in.SetPositive("dz", dz, err);
  in.SetPositive("dt", dt, err);
  if (dz != 0) dt_dz = dt / dz;
}
