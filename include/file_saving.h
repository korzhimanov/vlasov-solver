/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file file_saving.h
 * \brief The header file which defines some methods for convenient work with files
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#ifndef FILE_SAVING_H
#define FILE_SAVING_H

#include <cstdlib>
#include <cstdio>
#include <cstdarg>
#include <cstring>
#include <iostream>
#include <fstream>
#include <zlib.h>
#include <sys/stat.h>

namespace filesaving
{

FILE *open_file(const char *mode, const char *main_dir, const char *name, ...);
void close_file(FILE *file);

void create_dir(char *out, const char *main_dir, const char *dir, ...);
void create_dir(const char *main_dir, const char *dir, ...);

void save_file_1D(float *a, const int num, const char *main_dir, const char *name, ...);
void save_file_1D(float *a, const int num, const int step, const char *main_dir, const char *name, ...);
void save_file_1D(double *a, const int num, const char *main_dir, const char *name, ...);
void save_file_1D(double *a, const int num, const int step, const char *main_dir, const char *name, ...);
void save_file_1D(int *a, const int num, const char *main_dir, const char *name, ...);
void save_file_1D(int *a, const int num, const int step, const char *main_dir, const char *name, ...);

template<typename Type> void save_file_1D_bin(Type *a, const int num, const char *main_dir, const char *name, ...);
template<typename Type> void save_file_1D_gzip(Type *a, const int num, const char *main_dir, const char *name, ...);

void save_file_2D(double *a, int columns, int strings, const char *main_dir, const char *name, ...);
void save_file_2D_transpose(double *a, int columns, int strings, const char *main_dir, const char *name, ...);

} // namespace filesaving

template<typename Type>
void filesaving::save_file_1D_bin(Type *a, const int num, const char *main_dir, const char *name, ...)
{
    char temp[512];
    strcpy(temp, main_dir);
    va_list var; va_start(var, name);
    char name1[512];
    vsprintf(name1, name, var);
    strcat(temp, name1);

    std::ofstream out(temp, std::ios_base::binary|std::ios_base::app|std::ios_base::out);
    out.write((char*)a, num*sizeof(Type));
    out.close();

    va_end(var);
}

template<typename Type>
void filesaving::save_file_1D_gzip(Type *a, const int num, const char *main_dir, const char *name, ...)
{
    char temp[512];
    strcpy(temp, main_dir);
    va_list var; va_start(var, name);
    char name1[512];
    vsprintf(name1, name, var);
    strcat(temp, name1);

    gzFile out = gzopen(temp, "ab");
    gzwrite(out, (char*)a, num*sizeof(Type));
    gzclose(out);

    va_end(var);
}

#endif // FILE_SAVING_H
