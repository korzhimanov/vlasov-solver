/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file file_saving.h
 * \brief The header file which defines some methods for convenient work with
 * files
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#ifndef FILE_SAVING_H
#define FILE_SAVING_H

#include <cstdlib>
#include <fstream>
#include <zlib.h>

static std::ofstream fs;
static gzFile gz;

namespace filesaving {
void create_dir(std::string);

template <typename Type>
void save_file_1D_txt(Type *a, const int num, const std::string name);
template <typename Type>
void save_file_1D_txt(Type *a, const int num, const int step,
                      const std::string name);
template <typename Type>
void save_file_1D_bin(Type *a, const int num, const std::string name);
template <typename Type>
void save_file_1D_gzip(Type *a, const int num, const std::string name);

template <typename Type>
void save_file_2D_txt(Type *a, const int columns, const int rows,
                      const std::string name,
                      const bool need_transpose = false);
}  // namespace filesaving

template <typename Type>
void filesaving::save_file_1D_txt(Type *a, const int num,
                                  const std::string name) {
  fs.open(name.c_str(), std::ios_base::app | std::ios_base::out);

  int count;
  Type *b = a;
  for (count = 0; count < num; count++) fs << *(b++) << "\t";
  fs << std::endl;

  fs.close();
}

template <typename Type>
void filesaving::save_file_1D_txt(Type *a, const int num, const int step,
                                  const std::string name) {
  fs.open(name.c_str(), std::ios_base::app | std::ios_base::out);

  int count;
  for (count = 0; count < num * step; count += step) fs << *(a + count) << "\t";
  fs << std::endl;

  fs.close();
}

template <typename Type>
void filesaving::save_file_1D_bin(Type *a, const int num,
                                  const std::string name) {
  fs.open(name.c_str(),
          std::ios_base::binary | std::ios_base::app | std::ios_base::out);
  fs.write((char *)a, num * sizeof(Type));
  fs.close();
}

template <typename Type>
void filesaving::save_file_1D_gzip(Type *a, const int num,
                                   const std::string name) {
  gz = gzopen(name.c_str(), "ab");
  gzwrite(gz, (char *)a, num * sizeof(Type));
  gzclose(gz);
}

template <typename Type>
void filesaving::save_file_2D_txt(Type *a, const int columns, const int rows,
                                  const std::string name,
                                  const bool need_transpose = false) {
  fs.open(name.c_str(), std::ios_base::app | std::ios_base::out);

  int i, j, nx, ny;
  if (need_transpose) {
    nx = rows;
    ny = columns;
  } else {
    nx = columns;
    ny = rows;
  }
  for (i = 0; i < nx; i++) {
    for (j = 0; j < ny; j++) fs << *(a + j * nx + i) << "\t";
    fs << std::endl;
  }

  fs.close();
}

#endif  // FILE_SAVING_H
