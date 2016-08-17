/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file mesh.h
 * \brief The header file which defines Mesh class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#ifndef MESH_H
#define MESH_H

/**
 * /class Mesh
 * The class implements mesh parameters
*/
class Mesh
{
    public:
        double dz, dt, dt_dz;
        int MAX_Z, MAX_T;

    public:
        Mesh();
        virtual ~Mesh();
};

#endif // MESH_H
