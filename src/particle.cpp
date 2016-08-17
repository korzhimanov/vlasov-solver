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

#include "particle.h"
#include <cmath>

Particle::Particle()
{
    r = new double[3];
    p = new double[3];
    m = 1.;
    q = 1;
}

Particle::Particle(double* position, double* momentum, double mass, double charge)
{
    r = new double[3];
    p = new double[3];
    for (int i = 0; i < 3; i++)
    {
        r[i] = position[i];
        p[i] = momentum[i];
    }
    m = mass;
    q = int(charge);
}

Particle::~Particle()
{
    delete[] r;
    delete[] p;
}

void Particle::SetPosition(double *position)
{
    for (int i = 0; i < 3; i++)
        r[i] = position[i];
}

void Particle::SetMomentum(double *momentum)
{
    for (int i = 0; i < 3; i++)
        p[i] = momentum[i];
}

double Particle::GetGamma()
{
    return sqrt(1. + p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
}

void Particle::MakeStep(double dt, double *E, double *B)
{
    double F;
    double gamma = GetGamma();
    for (int i = 0; i < 3; i++)
    {
        F = E[i];
        /**
         * \todo Can be done smarter (and faster) without switch-case.
         * \todo Change division by gamma to multiplication to (1/gamma).
         */
        switch (i)
        {
            case 0 : F += (p[1]*B[2] - p[2]*B[1]) / gamma; break;
            case 1 : F += (p[2]*B[0] - p[0]*B[2]) / gamma; break;
            case 2 : F += (p[0]*B[1] - p[1]*B[0]) / gamma; break;
        }
        F *= -q/m;
        r[i] += (p[i] / (gamma*m) + .5 * F * dt) * dt;
        p[i] += F * dt;
    }
}
