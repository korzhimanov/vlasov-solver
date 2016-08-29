/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file particle.h
 * \brief The header file which defines Particle class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 * \todo Change Charge to double
 */

#ifndef INCLUDE_PARTICLE_H_
#define INCLUDE_PARTICLE_H_

/**
 * /class Particle
 * The class implements parameters of test particles and some useful methods on
 * them
*/
class Particle {
 public:
  Particle();
  Particle(double *position, double *momentum, double mass, double charge);
  ~Particle();
  void SetPosition(double *position);
  void SetMomentum(double *momentum);
  double GetGamma();
  void MakeStep(double dt, double *E, double *B);

  double *r;  // position in space
  double *p;  // momentum
  int q;      // charge
  double m;   // mass
};

#endif  // INCLUDE_PARTICLE_H_
