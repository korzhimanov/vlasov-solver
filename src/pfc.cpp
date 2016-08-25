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

#include "pfc.h"
#include <sstream>
#include <fstream>

PFC::PFC() : onesixth(1./6.)
{

}

PFC::PFC(int particle_type, pyinput *in, Mesh *m, double *n0, double *p0, int *err) : onesixth(1./6.)
{
    *err = Init(particle_type, in, m, n0, p0);
}

PFC::~PFC()
{
    delete profile;
    delete[] f1;
    delete[] f2;
    delete[] p;
    delete[] p2;
}

int PFC::Init(int particle_type, pyinput *in, Mesh *m, double *n0, double *p0)
{
    type = particle_type;
    mesh = m;
    N0 = *n0;
    P0 = *p0;

    std::stringstream ss;
    ss << type;
    std::string type_string(ss.str());

    profile = new pFunc(in->GetFunc("PROFILE_"+type_string));

    MAX_P = 64; if ( !in->SetPositive("MAX_P_"+type_string, &MAX_P) ) return 220;
    MASS = 1.; if ( !in->SetPositive("MASS_"+type_string, &MASS) ) return 230;
    CHARGE = 1.; if ( !in->Set("CHARGE_"+type_string, &CHARGE) ) return 240;
    q_m = CHARGE / MASS;
    q2_m2 = CHARGE * CHARGE / (MASS * MASS);
    MEAN_P = 0.; if ( !in->Set("MEAN_P_"+type_string, &MEAN_P) ) return 250;
    dp = 0.01; if ( !in->SetPositive("dp_"+type_string, &dp) ) return 260;

    if (fabs(MEAN_P) > 0.5*MAX_P*dp)
        std::cout << "WARNING! Mean momentum of " << type << " particle type is greater than maximal momentum of grid." << std::endl;
    if ((MEAN_P != 0.) && (fabs(MEAN_P) < dp))
        std::cout << "WARNING! Mean momentum of " << type << " particle type is less than momentum grid step." << std::endl;

    T_init =   0.02; if ( !in->SetPositive("T_init_"+type_string, &T_init) ) return 270;
    if (T_init < 0.5*dp*dp/MASS)
        std::cout << "WARNING! Temperature of " << type << " particle type cannot be resolved by momentum grid step." << std::endl;

    // some constants to speed up calculations
    halfq_m = .5 * CHARGE / MASS;
    half_qn0dz = .5  * CHARGE * N0 * mesh->dz * dp;
    quart_q2n0dtdp_m = .25 * CHARGE * CHARGE * N0 * mesh->dt * dp / MASS;
    qdt_mdp = CHARGE * mesh->dt / (MASS * dp);
    q2dt_dzdp = CHARGE * CHARGE * mesh->dt / (mesh->dz * dp);
    twicepisqrt3q2n0dpdzr_l = M_SQRT3 * CHARGE * CHARGE * N0 * dp * mesh->dz * M_PI * GSL_CONST_MKSA_ELECTRON_CHARGE * GSL_CONST_MKSA_ELECTRON_CHARGE / (GSL_CONST_MKSA_MASS_ELECTRON * GSL_CONST_MKSA_SPEED_OF_LIGHT * GSL_CONST_MKSA_SPEED_OF_LIGHT);

    return 0;
}

void PFC::AllocMemory()
{
    f1 = new double[mesh->MAX_Z*MAX_P];
    f2 = new double[mesh->MAX_Z*MAX_P];

    p  = new double[MAX_P];
    p2 = new double[MAX_P];
    for (int j = 0; j < MAX_P; j++)
    {
        p[j] = (j-MAX_P/2)*dp;
        p2[j] = p[j]*p[j];
    }
}

void PFC::SaveInput(std::ofstream& fs)
{
    fs << "\t\tCharge = " << CHARGE << " electron charge\n";
    fs << "\t\tMass = " << MASS << " electron mass\n";
    fs << "\t\tInitial temperature = " << T_init << " of the rest energy = " << T_init*511*MASS << " keV\n";
    fs << "\t\tMomentum cell size = " << dp << " Mc = " << dp*MASS*511 << " keV/c\n";
    fs << "\t\tTotal momentum cells = " << MAX_P << "\n";
    fs << "\t\tMaximal momentum = " << 0.5*MAX_P*dp << " Mc = " << 0.5*MAX_P*dp*MASS*511 << " keV/c\n";
}

void PFC::SetDistribution()
{
    // initialize plasma distribution
    maxf = 0.;
    double
        preexp = 1./( fabs(CHARGE) * sqrt(2.*M_PI*T_init) ),
        exponent = .5/T_init;
    for (int i = 0; i < mesh->MAX_Z; i++)
    {
        double A = profile->call(i) * preexp;
        for (int j = 0; j < MAX_P; j++)
            *(f1+i*MAX_P+j) = A * exp( -((j-MAX_P/2)*dp - MEAN_P) * ((j-MAX_P/2)*dp - MEAN_P) * exponent );
        if (maxf < *(f1+i*MAX_P+MAX_P/2))
            maxf = *(f1+i*MAX_P+MAX_P/2);
    }
}

void PFC::CalcLongitudinalField(double *ez)
{
    double *n = new double[mesh->MAX_Z];
    int i = 0,
        ii = 0,
        j = 0;
    #pragma omp parallel for private(i, ii, j)
    for (i = 0; i < mesh->MAX_Z; i++)
    {
        ii = i*MAX_P;
        n[i] = 0.;
        for (j = 0; j < MAX_P; j++)
            n[i] += *(f1+ii+j);
    }

    ez[0] -= half_qn0dz * n[0];
    for (i = 1; i < mesh->MAX_Z; i++)
        ez[i] -= half_qn0dz * (n[i] + n[i-1]);
    ez[mesh->MAX_Z] -= half_qn0dz * n[mesh->MAX_Z-1];

    delete[] n;
}

void PFC::CalcCurrent(FDTD *fdtd, double *ax, double *ay, double *a2)
{
    double current, a2i1;
    int i, ii;
    #pragma omp parallel for private(i, ii, current, a2i1)
    for (i = 1; i < mesh->MAX_Z; i++)
    {
        ii = i*MAX_P;
        a2i1 = 1. + ax[i]*ax[i]*q2_m2 + (ay[i]*q_m + P0)*(ay[i]*q_m + P0);
        current = 0.;
        for (int j = 0; j < MAX_P; j++)
        {
            current += (*(f1+ii+j) + *(f1+ii-MAX_P+j)) / sqrt(a2i1 + p2[j]);
        }
        current *= quart_q2n0dtdp_m;

        fdtd->ex[i] += (ax[i+1] + ax[i-1]) * current;
        fdtd->ey[i] += (ay[i+1] + ay[i-1] + 2*P0/q_m) * current;
    }
}

void PFC::MakeStep(double *ez, double *ax, double *ay)
{
    memset(f2, 0, sizeof(double)*mesh->MAX_Z*MAX_P);

    // calculate the distribution function at the first half-step
    int i, ii, j, new_cell;
    double a2i, a2i1, vdt_dz, flux_part, fold;
    #pragma omp parallel private(i, ii, a2i, a2i1, j, vdt_dz, new_cell, flux_part, fold)
    {
    #pragma omp for
    for (i = 0; i < mesh->MAX_Z; i++)
    {
        ii = i*MAX_P;
        a2i  = q2_m2 * ax[i]*ax[i] + (q_m*ay[i]+P0)*(q_m*ay[i]+P0);
        a2i1 = q2_m2 * ax[i+1]*ax[i+1] + (q_m*ay[i+1]+P0)*(q_m*ay[i+1]+P0);
        for (j = 0; j < MAX_P; j++)
        {
            if (*(f1+ii+j) > 1.e-14)
            {
                vdt_dz = - qdt_mdp * ez[i] - q2dt_dzdp * (sqrt(1. + a2i1 + p2[j]) - sqrt(1. + a2i + p2[j])); // calculate the relative velocity in phase space

                new_cell = j;
                flux_part = 0.;

                if (vdt_dz > 0.)
                {
                    new_cell = j - 1 + (int)ceil(vdt_dz);
                    flux_part = vdt_dz - floor(vdt_dz);
                    if (new_cell >= MAX_P - 1)
                    {
                        new_cell = MAX_P - 1;
                        flux_part = 0.;
                    }
                }
                else
                {
                    new_cell = j + (int)ceil(vdt_dz);
                    flux_part = vdt_dz - ceil(vdt_dz);
                    if (new_cell <= 0)
                    {
                        new_cell = 0;
                        flux_part = 0.;
                    }
                }

                fold = *(f1+ii+j);

                if (vdt_dz > 0.)
                {
                    if (j != MAX_P - 1)
                        fold += SlopeCorrectorRight(*(f1+ii+j), *(f1+ii+j+1)) * onesixth * (1.-flux_part) * (2.-flux_part) * (*(f1+ii+j+1) - *(f1+ii+j));
                    if (j != 0)
                        fold += SlopeCorrectorLeft (*(f1+ii+j), *(f1+ii+j-1)) * onesixth * (1.-flux_part) * (1.+flux_part) * (*(f1+ii+j) - *(f1+ii+j-1));

                    *(f2+ii+new_cell) += *(f1+ii+j) - flux_part*fold;

                    if (new_cell != MAX_P - 1)
                        *(f2+ii+new_cell+1) += flux_part*fold;
                }
                else
                {
                    if (j != MAX_P - 1)
                        fold -= SlopeCorrectorRight(*(f1+ii+j), *(f1+ii+j+1)) * onesixth * (1.-flux_part) * (1.+flux_part) * (*(f1+ii+j+1) - *(f1+ii+j));
                    if (j != 0)
                        fold -= SlopeCorrectorLeft (*(f1+ii+j), *(f1+ii+j-1)) * onesixth * (1.+flux_part) * (2.+flux_part) * (*(f1+ii+j) - *(f1+ii+j-1));

                    *(f2+ii+new_cell) += *(f1+ii+j) + flux_part*fold;

                    if (new_cell != 0)
                        *(f2+ii+new_cell-1) -= flux_part*fold;
                }
            }
        }
    }

    #pragma omp barrier
    #pragma omp single
    {
    memset(f1, 0, sizeof(double)*mesh->MAX_Z*MAX_P);
    }

// calculate the distribution function at the second half-step
    #pragma omp for
    for (i = 0; i < mesh->MAX_Z; i++)
    {
        ii = i*MAX_P;
        a2i  = q2_m2 * ax[i]*ax[i] + (q_m*ay[i]+P0)*(q_m*ay[i]+P0);
        a2i1 = q2_m2 * ax[i+1]*ax[i+1] + (q_m*ay[i+1]+P0)*(q_m*ay[i+1]+P0);
        for (j = 0; j < MAX_P; j++)
        {
            if (*(f2+ii+j) > 1.e-14)
            {
                vdt_dz = 2.*mesh->dt_dz*p[j]/(sqrt(1. + a2i1 + p2[j]) + sqrt(1. + a2i + p2[j])); // calculate the relative velocity in phase space

                new_cell = i;
                flux_part = 0.;

                if (vdt_dz > 0.)
                {
                    new_cell = i - 1 + (int)ceil(vdt_dz);
                    flux_part = vdt_dz - floor(vdt_dz);
                    if (new_cell >= mesh->MAX_Z - 1)
                    {
                        new_cell = mesh->MAX_Z - 1;
                        flux_part = 0.;
                    }
                }
                else
                {
                    new_cell = i + (int)ceil(vdt_dz);
                    flux_part = vdt_dz - ceil(vdt_dz);
                    if (new_cell <= 0)
                    {
                        new_cell = 0;
                        flux_part = 0.;
                    }
                }

                fold = *(f2+ii+j);
                if (vdt_dz > 0.)
                {
                    if (i != mesh->MAX_Z - 1)
                        fold += SlopeCorrectorRight(*(f2+ii+j), *(f2+ii+j+MAX_P)) * onesixth * (1.-flux_part) * (2.-flux_part) * (*(f2+ii+j+MAX_P) - *(f2+ii+j));
                    if (i != 0)
                        fold += SlopeCorrectorLeft (*(f2+ii+j), *(f2+ii+j-MAX_P)) * onesixth * (1.-flux_part) * (1.+flux_part) * (*(f2+ii+j) - *(f2+ii+j-MAX_P));

                    *(f1+new_cell*MAX_P+j) += *(f2+ii+j) - flux_part*fold;

                    if (new_cell != mesh->MAX_Z - 1)
                        *(f1+(new_cell+1)*MAX_P+j) += flux_part*fold;
                }
                else
                {
                    if (i != mesh->MAX_Z - 1)
                        fold -= SlopeCorrectorRight(*(f2+ii+j), *(f2+ii+j+MAX_P)) * onesixth * (1.-flux_part) * (1.+flux_part) * (*(f2+ii+j+MAX_P) - *(f2+ii+j));
                    if (i != 0)
                        fold -= SlopeCorrectorLeft(*(f2+ii+j), *(f2+ii+j-MAX_P)) * onesixth * (1.+flux_part) * (2.+flux_part) * (*(f2+ii+j) - *(f2+ii+j-MAX_P));

                    *(f1+new_cell*MAX_P+j) += *(f2+ii+j) + flux_part*fold;

                    if (new_cell != 0)
                        *(f1+(new_cell-1)*MAX_P+j) -= flux_part*fold;
                }
            }
        }
    }
    }
}

double PFC::KineticEnergy(double* ax, double *ay)
{
    return KineticEnergy(ax, ay, 0, mesh->MAX_Z-1);
}

double PFC::KineticEnergy(double* ax, double* ay, int left_bound, int right_bound)
{
    int i1, i2;
    if (left_bound < 0) i1 = 0;
    else i1 = left_bound;

    if (right_bound > mesh->MAX_Z-1) i2 = mesh->MAX_Z-1;
    else i2 = right_bound;

    if (i2 < i1) return 0.;

    double energy = 0.;
    for (int i = i1; i < i2; i++)
    {
        int ii = i*MAX_P;
        double a2i  = q2_m2 * ax[i]*ax[i] + (q_m*ay[i]+P0)*(q_m*ay[i]+P0);
        double a2i1 = q2_m2 * ax[i+1]*ax[i+1] + (q_m*ay[i+1]+P0)*(q_m*ay[i+1]+P0);
        for (int j = 0; j < MAX_P; j++)
            energy += (sqrt(1. + a2i + p2[j]) + sqrt(1. + a2i1 + p2[j]) - 2.) * *(f1+ii+j);
    }

    return energy * .5 * fabs(CHARGE) * MASS * N0 * dp * mesh->dz;
}

void PFC::SaveConcentrationTxt(std::string filename)
{
    float *n = new float[mesh->MAX_Z];

    CalcConcDstr(n);

    filesaving::save_file_1D_txt(n, mesh->MAX_Z, filename+".txt");

    delete[] n;
}

void PFC::SaveConcentrationBin(std::string filename)
{
    float *n = new float[mesh->MAX_Z];

    CalcConcDstr(n);

    std::ofstream out((filename+".bin").c_str(), std::ios_base::binary|std::ios_base::app);
    out.write((char*)n, mesh->MAX_Z*sizeof(float));
    out.close();

    delete[] n;
}

void PFC::SaveConcentrationGZip(std::string filename)
{
    float *n = new float[mesh->MAX_Z];

    CalcConcDstr(n);

    gzFile out = gzopen((filename+".gz").c_str(), "ab");
    gzwrite(out, (char*)n, mesh->MAX_Z*sizeof(float));
    gzclose(out);

    delete[] n;
}

void PFC::SaveDstrFunctionTxt(std::string filename)
{
    filesaving::save_file_2D_txt(f1, MAX_P, mesh->MAX_Z, filename+".txt", true);
}

void PFC::SaveDstrFunctionBin(std::string filename)
{
    std::ofstream out((filename+".bin").c_str(), std::ios_base::binary);
    out.write((char*)f1, MAX_P*mesh->MAX_Z*sizeof(double));
    out.close();
}

void PFC::SaveDstrFunctionGZip(std::string filename)
{
    gzFile out = gzopen((filename+".gz").c_str(), "wb");
    gzwrite(out, (char*)f1, MAX_P*mesh->MAX_Z*sizeof(double));
    gzclose(out);
}

void PFC::SaveEmitRadiationTxt(std::string filename, int max_harmonic, int dn, double dt, FDTD *fdtd, double *ez, double *ax, double *ay, double *a2)
{
    double *I = new double[max_harmonic/dn+1];
    I[0] = 0.;

    for (int n = dn; n <= max_harmonic; n+=dn)
    {
        I[n/dn] = CalcEmitRadiation(double(n), fdtd, ez, ax, ay, a2)*dt;
    }

    filesaving::save_file_1D_txt(I, max_harmonic/dn+1, filename+".txt");

    delete[] I;
}

void PFC::SaveEmitRadiationBin(std::string filename, int max_harmonic, int dn, double dt, FDTD *fdtd, double *ez, double *ax, double *ay, double *a2)
{
    double *I = new double[max_harmonic/dn+1];
    I[0] = 0.;

    for (int n = dn; n <= max_harmonic; n+=dn)
    {
        I[n/dn] = CalcEmitRadiation(double(n), fdtd, ez, ax, ay, a2)*dt;
    }

    std::ofstream out((filename+".bin").c_str(), std::ios_base::binary);
    out.write((char*)I, (max_harmonic/dn+1)*sizeof(double));
    out.close();

    delete[] I;
}

void PFC::CalcConcDstr(float* n)
{
    memset(n, 0, sizeof(float)*mesh->MAX_Z);
    for (int i = 0; i < mesh->MAX_Z; i++)
    {
        int ii = i*MAX_P;
        for (int j = 0; j < MAX_P; j++)
            n[i] += float(*(f1+ii+j));
        n[i] *= float(N0*dp);
    }
}

double PFC::CalcEmitRadiation(double freq, FDTD *fdtd, double *ez, double *ax, double *ay, double *a2)
{
    double ae[3]; // acceleration by electric field
    double ab[3]; // acceleration by magnetic field
    double an[3]; // normal acceleration
    double v[3]; // velocity

    double gamma; // gamma-factor
    double gamma2; // gamma-factor squared
    double inv_gamma; // inverse gamma-factor
    double av; // scalar product of acceleration and velocity
    double gamma2_curv_rad; // gamma squared divided by curvature radius
    double cut_freq; // cutoff frequency
    double dI_dt = 0.; // power of radiation per unit of frequency
    for (int i = 0; i < mesh->MAX_Z; i++)
    {
        ae[0] = -halfq_m * (fdtd->ex[i]+fdtd->ex[i+1]);
        ae[1] = -halfq_m * (fdtd->ey[i]+fdtd->ey[i+1]);
        ae[2] = -halfq_m * (ez[i]+ez[i+1]);

        int ii = i*MAX_P;
        double a2_1 = 1. + a2[i]*q2_m2;
        double px = halfq_m * (ax[i] + ax[i+1]);
        double py = halfq_m * (ay[i] + ay[i+1]);
        for (int j = 0; j < MAX_P; j++)
        {
            if (*(f1+ii+j) > 0.)
            {
                gamma2 = a2_1 + p2[j];
                if (gamma2 > 2.)
                {
                    gamma = sqrt(gamma2);
                    inv_gamma = 1./gamma;

                    v[0] = px   * inv_gamma;
                    v[1] = py   * inv_gamma;
                    v[2] = p[j] * inv_gamma;

                    ab[0] =  q_m *  v[2] * fdtd->hy[i];
                    ab[1] = -q_m *  v[2] * fdtd->hx[i];
                    ab[2] = -q_m * (v[0] * fdtd->hy[i] - v[1] * fdtd->hx[i]);

                    av = ae[0]*v[0] + ae[1]*v[1] + ae[2]*v[2];

                    an[0] = inv_gamma * (ae[0] + ab[0] - av*v[0]);
                    an[1] = inv_gamma * (ae[1] + ab[1] - av*v[1]);
                    an[2] = inv_gamma * (ae[2] + ab[2] - av*v[2]);

                    gamma2_curv_rad = sqrt(an[0]*an[0]+an[1]*an[1]+an[2]*an[2])*gamma2*gamma2/(gamma2-1.);

                    cut_freq = 3.*gamma*gamma2_curv_rad;

                    if (cut_freq*mesh->dt > gamma2 && freq < 10.*cut_freq)
                    {
                        dI_dt += gamma2_curv_rad * gsl_sf_synchrotron_1(2.*freq/cut_freq) * *(f1+ii+j);
                    }
                }
            }
        }
    }

    return twicepisqrt3q2n0dpdzr_l*dI_dt;
}
