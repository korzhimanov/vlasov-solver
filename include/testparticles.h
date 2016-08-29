/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file testparticles.h
 * \brief The header file which defines TestParticles class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#ifndef INCLUDE_TESTPARTICLES_H_
#define INCLUDE_TESTPARTICLES_H_

#include "include/mesh.h"
#include "include/particle.h"
#include "include/pyinput.h"

class TestParticles {
 public:
  Particle *prt;  // particles
  int particles_number, start_point;
  double interval, mass, charge,
      mean_initial_momentum;  // initial mean momentum along y-axis (useful for
                              // simulation of oblique irradiation in boosted
                              // frame)

 public:
  TestParticles(PyInput &, const double &theta, int &err);
  virtual ~TestParticles();
  void AllocMemory();
  void InitParticles(Mesh *);
};

#endif  // INCLUDE_TESTPARTICLES_H_
