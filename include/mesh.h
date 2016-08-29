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

#include "pyinput.h"

/**
 * /class Mesh
 * The class to store and initialeze mesh parameters
*/
class Mesh {
 public:
  int ppw, MAX_Z, MAX_T;
  double dz, dt, dt_dz;

 public:
  Mesh();
  Mesh(pyinput*, int*);
  int Init(pyinput*);
};

#endif  // MESH_H
