/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file mesh.h
 * \brief The header file which defines Mesh struct
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#ifndef MESH_H
#define MESH_H

/**
 * /struct Mesh
 * The structure for mesh parameters
*/
struct Mesh
{
    double dz, dt, dt_dz;
    int MAX_Z, MAX_T;
};

#endif // MESH_H
