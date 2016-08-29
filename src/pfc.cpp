/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file pfc.cpp
 * \brief The source file which implements the methods of the PFC class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#include "include/pfc.h"

#include <zlib.h>

#include <iostream>
#include <sstream>

#include "include/file_saving.h"
#include "include/mymath.h"

PFC::PFC() {}

PFC::PFC(int particle_type, pyinput *in, Mesh *m, double *n0, double *p0,
         int *err) {
  *err = Init(particle_type, in, m, n0, p0);
}

PFC::~PFC() {
  delete profile;
  delete[] f1;
  delete[] f2;
  delete[] conc;
  delete[] p;
  delete[] p2;
}

int PFC::Init(int particle_type, pyinput *in, Mesh *m, double *n0, double *p0) {
  type = particle_type;
  mesh = m;
  N0 = *n0;
  P0 = *p0;

  std::stringstream ss;
  ss << type;
  std::string type_string(ss.str());

  profile = new pFunc(in->GetFunc("PROFILE_" + type_string));

  MAX_P = 64;
  if (!in->SetPositive("MAX_P_" + type_string, &MAX_P))
    return 220;
  MASS = 1.;
  if (!in->SetPositive("MASS_" + type_string, &MASS))
    return 230;
  CHARGE = 1.;
  if (!in->Set("CHARGE_" + type_string, &CHARGE))
    return 240;
  q_m = CHARGE / MASS;
  q2_m2 = CHARGE * CHARGE / (MASS * MASS);
  MEAN_P = 0.;
  if (!in->Set("MEAN_P_" + type_string, &MEAN_P))
    return 250;
  dp = 0.01;
  if (!in->SetPositive("dp_" + type_string, &dp))
    return 260;

  if (fabs(MEAN_P) > 0.5 * MAX_P * dp)
    std::cout << "WARNING! Mean momentum of " << type
              << " particle type is greater than maximal momentum of grid."
              << std::endl;
  if ((MEAN_P != 0.) && (fabs(MEAN_P) < dp))
    std::cout << "WARNING! Mean momentum of " << type
              << " particle type is less than momentum grid step." << std::endl;

  T_init = 0.02;
  if (!in->SetPositive("T_init_" + type_string, &T_init))
    return 270;
  if (T_init < 0.5 * dp * dp / MASS)
    std::cout << "WARNING! Temperature of " << type
              << " particle type cannot be resolved by momentum grid step."
              << std::endl;

  // some constants to speed up calculations
  n0dp = N0 * dp;
  halfq_m = .5 * CHARGE / MASS;
  halfqdz = .5 * CHARGE * mesh->dz;
  quart_q2n0dtdp_m = .25 * CHARGE * CHARGE * N0 * mesh->dt * dp / MASS;
  qdt_mdp = CHARGE * mesh->dt / (MASS * dp);
  q2dt_dzdp = CHARGE * CHARGE * mesh->dt / (mesh->dz * dp);

  return 0;
}

void PFC::AllocMemory() {
  f1 = new double[mesh->MAX_Z * MAX_P];
  f2 = new double[mesh->MAX_Z * MAX_P];

  conc = new double[mesh->MAX_Z];

  p = new double[MAX_P];
  p2 = new double[MAX_P];
  for (int j = 0; j < MAX_P; j++) {
    p[j] = (j - MAX_P / 2) * dp;
    p2[j] = p[j] * p[j];
  }
}

void PFC::SaveInput(std::ofstream &fs) {
  fs << "\t\tCharge = " << CHARGE << " electron charge\n";
  fs << "\t\tMass = " << MASS << " electron mass\n";
  fs << "\t\tInitial temperature = " << T_init
     << " of the rest energy = " << T_init * 511 * MASS << " keV\n";
  fs << "\t\tMomentum cell size = " << dp << " Mc = " << dp * MASS * 511
     << " keV/c\n";
  fs << "\t\tTotal momentum cells = " << MAX_P << "\n";
  fs << "\t\tMaximal momentum = " << 0.5 * MAX_P * dp
     << " Mc = " << 0.5 * MAX_P * dp * MASS * 511 << " keV/c\n";
}

void PFC::SetDistribution() {
  // initialize plasma distribution
  maxf = 0.;
  double preexp = 1. / (fabs(CHARGE) * sqrt(2. * M_PI * T_init)),
         exponent = .5 / T_init;
  for (int i = 0; i < mesh->MAX_Z; i++) {
    double A = profile->call(i) * preexp;
    for (int j = 0; j < MAX_P; j++)
      *(f1 + i * MAX_P + j) =
          A * exp(-((j - MAX_P / 2) * dp - MEAN_P) *
                  ((j - MAX_P / 2) * dp - MEAN_P) * exponent);
    if (maxf < *(f1 + i * MAX_P + MAX_P / 2))
      maxf = *(f1 + i * MAX_P + MAX_P / 2);
  }

  CalcConcDstr();
}

void PFC::CalcLongitudinalField(double *ez) {
  CalcConcDstr();

  ez[0] -= halfqdz * conc[0];
  for (int i = 1; i < mesh->MAX_Z; i++)
    ez[i] -= halfqdz * (conc[i] + conc[i - 1]);
  ez[mesh->MAX_Z] -= halfqdz * conc[mesh->MAX_Z - 1];
}

void PFC::CalcCurrent(FDTD *fdtd, double *ax, double *ay, double *a2) {
  double current, a2i1;
  int i, ii;
#pragma omp parallel for private(i, ii, current, a2i1)
  for (i = 1; i < mesh->MAX_Z; i++) {
    ii = i * MAX_P;
    a2i1 = 1. + ax[i] * ax[i] * q2_m2 + (ay[i] * q_m + P0) * (ay[i] * q_m + P0);
    current = 0.;
    for (int j = 0; j < MAX_P; j++) {
      current += (*(f1 + ii + j) + *(f1 + ii - MAX_P + j)) / sqrt(a2i1 + p2[j]);
    }
    current *= quart_q2n0dtdp_m;

    fdtd->ex[i] += (ax[i + 1] + ax[i - 1]) * current;
    fdtd->ey[i] += (ay[i + 1] + ay[i - 1] + 2 * P0 / q_m) * current;
  }
}

void PFC::MakeStep(double *ez, double *ax, double *ay) {
  memset(f2, 0, sizeof(double) * mesh->MAX_Z * MAX_P);

  // calculate the distribution function at the first half-step
  int i, ii, j, new_cell;
  double a2i, a2i1, vdt_dz, flux_part, fold;
#pragma omp parallel private(i, ii, a2i, a2i1, j, vdt_dz, new_cell, flux_part, \
                             fold)
  {
#pragma omp for
    for (i = 0; i < mesh->MAX_Z; i++)
      if (conc[i] > 1e-14) {
        ii = i * MAX_P;
        a2i = q2_m2 * ax[i] * ax[i] + (q_m * ay[i] + P0) * (q_m * ay[i] + P0);
        a2i1 = q2_m2 * ax[i + 1] * ax[i + 1] +
               (q_m * ay[i + 1] + P0) * (q_m * ay[i + 1] + P0);
        for (j = 0; j < MAX_P; j++) {
          if (*(f1 + ii + j) > 1.e-14) {
            // calculate the relative velocity in phase space
            vdt_dz =
                -qdt_mdp * ez[i] -
                q2dt_dzdp * (sqrt(1. + a2i1 + p2[j]) - sqrt(1. + a2i + p2[j]));

            new_cell = j;
            flux_part = 0.;

            if (vdt_dz > 0.) {
              new_cell = j - 1 + static_cast<int>(ceil(vdt_dz));
              flux_part = vdt_dz - floor(vdt_dz);
              if (new_cell >= MAX_P - 1) {
                new_cell = MAX_P - 1;
                flux_part = 0.;
              }
            } else {
              new_cell = j + static_cast<int>(ceil(vdt_dz));
              flux_part = vdt_dz - ceil(vdt_dz);
              if (new_cell <= 0) {
                new_cell = 0;
                flux_part = 0.;
              }
            }

            fold = *(f1 + ii + j);

            if (vdt_dz > 0.) {
              if (j != MAX_P - 1)
                fold += CorrectedRightSlope(f1 + ii + j, f1 + ii + j + 1) *
                        mymath::onesixth * (1. - flux_part) * (2. - flux_part);
              if (j != 0)
                fold += CorrectedLeftSlope(f1 + ii + j, f1 + ii + j - 1) *
                        mymath::onesixth * (1. - flux_part) * (1. + flux_part);

              fold *= flux_part;

              *(f2 + ii + new_cell) += *(f1 + ii + j) - fold;

              if (new_cell != MAX_P - 1)
                *(f2 + ii + new_cell + 1) += fold;
            } else {
              if (j != MAX_P - 1)
                fold -= CorrectedRightSlope(f1 + ii + j, f1 + ii + j + 1) *
                        mymath::onesixth * (1. - flux_part) * (1. + flux_part);
              if (j != 0)
                fold -= CorrectedLeftSlope(f1 + ii + j, f1 + ii + j - 1) *
                        mymath::onesixth * (1. + flux_part) * (2. + flux_part);

              fold *= flux_part;

              *(f2 + ii + new_cell) += *(f1 + ii + j) + fold;

              if (new_cell != 0)
                *(f2 + ii + new_cell - 1) -= fold;
            }
          }
        }
      }

#pragma omp barrier
#pragma omp single
    { memset(f1, 0, sizeof(double) * mesh->MAX_Z * MAX_P); }

// calculate the distribution function at the second half-step
#pragma omp for
    for (i = 0; i < mesh->MAX_Z; i++)
      if (conc[i] > 1e-14) {
        ii = i * MAX_P;
        a2i = q2_m2 * ax[i] * ax[i] + (q_m * ay[i] + P0) * (q_m * ay[i] + P0);
        a2i1 = q2_m2 * ax[i + 1] * ax[i + 1] +
               (q_m * ay[i + 1] + P0) * (q_m * ay[i + 1] + P0);
        for (j = 0; j < MAX_P; j++) {
          if (*(f2 + ii + j) > 1.e-14) {
            // calculate the relative velocity in phase space
            vdt_dz = 2. * mesh->dt_dz * p[j] /
                     (sqrt(1. + a2i1 + p2[j]) + sqrt(1. + a2i + p2[j]));

            new_cell = i;
            flux_part = 0.;

            if (vdt_dz > 0.) {
              new_cell = i - 1 + static_cast<int>(ceil(vdt_dz));
              flux_part = vdt_dz - floor(vdt_dz);
              if (new_cell >= mesh->MAX_Z - 1) {
                new_cell = mesh->MAX_Z - 1;
                flux_part = 0.;
              }
            } else {
              new_cell = i + static_cast<int>(ceil(vdt_dz));
              flux_part = vdt_dz - ceil(vdt_dz);
              if (new_cell <= 0) {
                new_cell = 0;
                flux_part = 0.;
              }
            }

            fold = *(f2 + ii + j);
            if (vdt_dz > 0.) {
              if (i != mesh->MAX_Z - 1)
                fold += CorrectedRightSlope(f2 + ii + j, f2 + ii + j + MAX_P) *
                        mymath::onesixth * (1. - flux_part) * (2. - flux_part);
              if (i != 0)
                fold += CorrectedLeftSlope(f2 + ii + j, f2 + ii + j - MAX_P) *
                        mymath::onesixth * (1. - flux_part) * (1. + flux_part);

              fold *= flux_part;

              *(f1 + new_cell * MAX_P + j) += *(f2 + ii + j) - fold;

              if (new_cell != mesh->MAX_Z - 1)
                *(f1 + (new_cell + 1) * MAX_P + j) += fold;
            } else {
              if (i != mesh->MAX_Z - 1)
                fold -= CorrectedRightSlope(f2 + ii + j, f2 + ii + j + MAX_P) *
                        mymath::onesixth * (1. - flux_part) * (1. + flux_part);
              if (i != 0)
                fold -= CorrectedLeftSlope(f2 + ii + j, f2 + ii + j - MAX_P) *
                        mymath::onesixth * (1. + flux_part) * (2. + flux_part);

              fold *= flux_part;

              *(f1 + new_cell * MAX_P + j) += *(f2 + ii + j) + fold;

              if (new_cell != 0)
                *(f1 + (new_cell - 1) * MAX_P + j) -= fold;
            }
          }
        }
      }
  }
}

// slope corrector for right boundary
double PFC::CorrectedRightSlope(double *centralvalue, double *rightvalue) {
  double slope = *rightvalue - *centralvalue, slope_c;
  if (slope < -1e-14) {
    slope_c = -2. * (maxf - *centralvalue);
    if (slope_c > slope)
      return slope_c;
    else
      return slope;
  }
  if (slope > 1e-14) {
    slope_c = 2. * (*centralvalue);
    if (slope_c < slope)
      return slope_c;
    else
      return slope;
  }
  return slope;
}

// slope corrector for left boundary
double PFC::CorrectedLeftSlope(double *centralvalue, double *leftvalue) {
  double slope = *centralvalue - *leftvalue, slope_c;
  if (slope > 1e-14) {
    slope_c = 2. * (maxf - *centralvalue);
    if (slope_c < slope)
      return slope_c;
    else
      return slope;
  }
  if (slope < -1e-14) {
    slope_c = -2. * (*centralvalue);
    if (slope_c > slope)
      return slope_c;
    else
      return slope;
  }
  return slope;
}

double PFC::KineticEnergy(double *ax, double *ay) {
  return KineticEnergy(ax, ay, 0, mesh->MAX_Z - 1);
}

double PFC::KineticEnergy(double *ax, double *ay, int left_bound,
                          int right_bound) {
  int i1, i2;
  if (left_bound < 0)
    i1 = 0;
  else
    i1 = left_bound;

  if (right_bound > mesh->MAX_Z - 1)
    i2 = mesh->MAX_Z - 1;
  else
    i2 = right_bound;

  if (i2 < i1)
    return 0.;

  double energy = 0.;
  for (int i = i1; i < i2; i++) {
    int ii = i * MAX_P;
    double a2i =
        q2_m2 * ax[i] * ax[i] + (q_m * ay[i] + P0) * (q_m * ay[i] + P0);
    double a2i1 = q2_m2 * ax[i + 1] * ax[i + 1] +
                  (q_m * ay[i + 1] + P0) * (q_m * ay[i + 1] + P0);
    for (int j = 0; j < MAX_P; j++)
      energy += (sqrt(1. + a2i + p2[j]) + sqrt(1. + a2i1 + p2[j]) - 2.) *
                *(f1 + ii + j);
  }

  return energy * .5 * fabs(CHARGE) * MASS * N0 * dp * mesh->dz;
}

void PFC::SaveConcentrationTxt(std::string filename) {
  filesaving::save_file_1D_txt(conc, mesh->MAX_Z, filename + ".txt");
}

void PFC::SaveConcentrationBin(std::string filename) {
  std::ofstream out((filename + ".bin").c_str(),
                    std::ios_base::binary | std::ios_base::app);
  out.write(reinterpret_cast<char *>(conc), mesh->MAX_Z * sizeof(double));
  out.close();
}

void PFC::SaveConcentrationGZip(std::string filename) {
  gzFile out = gzopen((filename + ".gz").c_str(), "ab");
  gzwrite(out, reinterpret_cast<char *>(conc), mesh->MAX_Z * sizeof(double));
  gzclose(out);
}

void PFC::SaveDstrFunctionTxt(std::string filename) {
  filesaving::save_file_2D_txt(f1, MAX_P, mesh->MAX_Z, filename + ".txt", true);
}

void PFC::SaveDstrFunctionBin(std::string filename) {
  std::ofstream out((filename + ".bin").c_str(), std::ios_base::binary);
  out.write(reinterpret_cast<char *>(f1), MAX_P * mesh->MAX_Z * sizeof(double));
  out.close();
}

void PFC::SaveDstrFunctionGZip(std::string filename) {
  gzFile out = gzopen((filename + ".gz").c_str(), "wb");
  gzwrite(out, reinterpret_cast<char *>(f1),
          MAX_P * mesh->MAX_Z * sizeof(double));
  gzclose(out);
}

void PFC::CalcConcDstr() {
  memset(conc, 0, sizeof(double) * mesh->MAX_Z);
  int i = 0, ii = 0, j = 0;
#pragma omp parallel for private(i, ii, j)
  for (i = 0; i < mesh->MAX_Z; i++) {
    ii = i * MAX_P;
    for (j = 0; j < MAX_P; j++)
      conc[i] += *(f1 + ii + j);
    conc[i] *= n0dp;
  }
}
