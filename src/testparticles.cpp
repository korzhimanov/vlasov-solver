/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file testparticles.cpp
 * \brief The source file which implements the methods of the TestParticles
 * class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#include "include/testparticles.h"

#include <iostream>

#include "include/errors.h"

TestParticles::TestParticles(pyinput &in, const double &theta, int &err)
    : particles_number(0),
      start_point(0),
      interval(1.),
      mass(1.),
      charge(0.),
      mean_initial_momentum(tan(theta)) {
  in.SetNotNegative("NUM_PRT", particles_number, err);

  if (particles_number > 0) {
    in.SetNotNegative("start_point", start_point, err);
    in.SetPositive("interval", interval, err);
    in.SetPositive("MASS_PRT", mass, err);
    in.Set("CHARGE_PRT", charge, err);
  }
}

TestParticles::~TestParticles() {
  if (particles_number > 0) delete[] prt;
}

void TestParticles::AllocMemory() {
  if (particles_number > 0) prt = new Particle[particles_number];
}

void TestParticles::InitParticles(Mesh *mesh) {
  if (particles_number > 0) {
    for (int i = 0; i < particles_number; i++) {
      prt[i].q = static_cast<int>(charge);
      prt[i].m = mass;
      for (int j = 0; j < 3; j++) {
        prt[i].r[j] = 0.;
        prt[i].p[j] = 0.;
      }
      prt[i].p[1] = mean_initial_momentum;
      prt[i].r[2] = (start_point + interval * i) * mesh->dz;
    }
  }
}
