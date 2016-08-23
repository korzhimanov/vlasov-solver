/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file file_saving.cpp
 * \brief The source file which implements some methods for convenient work with files
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#include "file_saving.h"
#include <string>
#include <sstream>

namespace filesaving {

FILE* open_file(const char *mode, const std::string main_dir, const char *name, ...)
{
    char temp[512];
    strcpy(temp, main_dir.c_str());
    va_list var; va_start(var, name);
    char name1[512];
    vsprintf(name1, name, var);
    strcat(temp, name1);

    FILE *fp;
    fp = fopen(temp, mode);
    if (!fp) std::cout << "ERROR: Cannot open file " << temp << std::endl;

    va_end(var);

    return fp;
}

void close_file(FILE *file)
{
    fclose(file);
}

// when creating directory put "/" at the end of its name
void create_dir(std::string dn)
{
    std::stringstream command;
#ifdef _WIN32
    command << "del ";
#else
    command << "rm -rf ";
#endif

    command << dn;
#ifdef _WIN32
    command << " /q";
#endif

    std::cout << command.str() << std::endl;

    system(command.str().c_str());

    command.str( std::string() );
    command.clear();
#ifdef _WIN32
    command << "mkdir " << dn;
    std::cout << command.str() << std::endl;
    system(command.str().c_str());
#else
    command << dn;
    std::cout << "mkdir " << command.str() << std::endl;
    mkdir(command.str().c_str(), 0755);
#endif
}

void save_file_1D(float *a, const int num, const std::string main_dir, const char *name, ...)
{
    va_list var; va_start(var, name);

    char temp[256];
    vsprintf(temp, name, var);

    FILE *file = open_file("a+", main_dir, temp);

    int count;
    float *b = a;
    for (count = 0; count < num; count++)
        fprintf(file, "%8.5f\t", *(b++));
    fprintf(file, "\n");

    close_file(file);

    va_end(var);
}

void save_file_1D(float *a, const int num, const int step, const std::string main_dir, const char *name, ...)
{
    va_list var; va_start(var, name);

    char temp[256];
    vsprintf(temp, name, var);

    FILE *file = open_file("a+", main_dir, temp);

    int count, until = num*step;
    for (count = 0; count < until; count += step)
        fprintf(file, "%8.5f\t", *(a+count));
    fprintf(file, "\n");

    close_file(file);

    va_end(var);
}

void save_file_1D(double *a, const int num, const std::string main_dir, const char *name, ...)
{
    va_list var; va_start(var, name);

    char temp[256];
    vsprintf(temp, name, var);

    FILE *file = open_file("a+", main_dir, temp);

    int count;
    double *b = a;
    for (count = 0; count < num; count++)
        fprintf(file, "%8.2e\t", *(b++));
    fprintf(file, "\n");

    close_file(file);

    va_end(var);
}

void save_file_1D(double *a, const int num, const int step, const std::string main_dir, const char *name, ...)
{
    va_list var; va_start(var, name);

    char temp[256];
    vsprintf(temp, name, var);

    FILE *file = open_file("a+", main_dir, temp);

    int count, until = num*step;
    for (count = 0; count < until; count += step) fprintf(file, "%8.2e\t", *(a+count));
    fprintf(file, "\n");

    close_file(file);

    va_end(var);
}

void save_file_1D(int *a, const int num, const std::string main_dir, const char *name, ...)
{
    va_list var; va_start(var, name);

    char temp[256];
    vsprintf(temp, name, var);

    FILE *file = open_file("a+", main_dir, temp);

    int count;
    int *b = a;
    for (count = 0; count < num; count++) fprintf(file, "%8d\t", *(b++));
    fprintf(file, "\n");

    close_file(file);

    va_end(var);
}

void save_file_1D(int *a, const int num, const int step, const std::string main_dir, const char *name, ...)
{
    va_list var; va_start(var, name);

    char temp[256];
    vsprintf(temp, name, var);

    FILE *file = open_file("a+",main_dir, temp);

    int count, until = num*step;
    for (count = 0; count < until; count += step) fprintf(file, "%8d\t", *(a+count));
    fprintf(file, "\n");

    close_file(file);

    va_end(var);
}

void save_file_2D(double *a, int columns, int strings, const std::string main_dir, const char *name, ...)
{
    va_list var; va_start(var, name);

    char temp[256];
    vsprintf(temp, name, var);

    FILE *file = open_file("a+",main_dir, temp);

    int count1, count2;
    double *b = a;
    for (count1 = 0; count1 < strings; count1++)
    {
        for (count2 = 0; count2 < columns; count2++)
            fprintf(file, "%8.5f\t", *(b++));
        fprintf(file, "\n");
    }

    close_file(file);

    va_end(var);
}

void save_file_2D_transpose(double *a, int columns, int strings, const std::string main_dir, const char *name, ...)
{
    va_list var; va_start(var, name);

    char temp[256];
    vsprintf(temp, name, var);

    FILE *file = open_file("a+", main_dir, temp);

    int count1, count2;
    for (count1 = 0; count1 < columns; count1++)
    {
        for (count2 = 0; count2 < strings; count2++)
            fprintf(file, "%8.5f\t", *(a + count2*columns + count1));
        fprintf(file, "\n");
    }

    close_file(file);

    va_end(var);
}

} // namespace filesaving
