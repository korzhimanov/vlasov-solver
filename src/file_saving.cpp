/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file file_saving.cpp
 * \brief The source file which implements the methods of the FileSaving class
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#include"file_saving.h"

FileSaving::FileSaving()
{
}

FileSaving::~FileSaving()
{
}

FILE* FileSaving::open_file(const char *mode, const char *main_dir, const char *name, ...)
{
    char temp[512];
    strcpy(temp, main_dir);
    va_list var; va_start(var, name);
    char name1[512];
    vsprintf(name1, name, var);
    strcat(temp, name1);

    FILE *fp;
    fp = fopen(temp, mode);
    if (!fp) cout << "ERROR: Cannot open file " << temp << endl;

    va_end(var);

    return fp;
}

void FileSaving::close_file(FILE *file)
{
    fclose(file);
}

// when creating directory put "/" at the end of its name
void FileSaving::create_dir(char *out, const char *main_dir, const char *dir, ...)
{
    char command[512];
#if SYSTEM==0
    sprintf(command, "del ");
#else
    sprintf(command, "rm -rf ");
#endif

    va_list var; va_start(var, dir);
    char temp[256];
    vsprintf(temp, dir, var);

    strcat(command, main_dir);
    strcat(command, temp);
#if SYSTEM==0
    {
        char args[10] = " /q";
        strcat(command, args);
    }
#endif

    cout << command << endl;

    system(command);

#if SYSTEM==0
    sprintf(command, "mkdir ");
    strcat(command, main_dir);
#else
    strcpy(command, main_dir);
#endif
    strcat(command, temp);

#if SYSTEM==0
    cout << command << endl;
#else
    cout << "mkdir " << command << endl;
#endif

#if SYSTEM==0
    system(command);
#else
    mkdir(command, 0777);
#endif

    strcpy(out, temp);

    va_end(var);
}

// when creating directory put "/" at the end of its name
void FileSaving::create_dir(const char *main_dir, const char *dir, ...)
{
    char command[512];
#if SYSTEM==0
    sprintf(command, "del ");
#else
    sprintf(command, "rm -rf ");
#endif

    va_list var; va_start(var, dir);
    char temp[256];
    vsprintf(temp, dir, var);

    strcat(command, main_dir);
    strcat(command, temp);
#if SYSTEM==0
    {
        char args[10] = " /q";
        strcat(command, args);
    }
#endif

    cout << command << endl;

    system(command);

#if SYSTEM==0
    sprintf(command, "mkdir ");
    strcat(command, main_dir);
#else
    strcpy(command, main_dir);
#endif
    strcat(command, temp);

#if SYSTEM==0
    cout << command << endl;
#else
    cout << "mkdir " << command << endl;
#endif

#if SYSTEM==0
    system(command);
#else
    mkdir(command, 0777);
#endif

    va_end(var);
}

void FileSaving::save_file_1D(float *a, const int num, const char *main_dir, const char *name, ...)
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

void FileSaving::save_file_1D(float *a, const int num, const int step, const char *main_dir, const char *name, ...)
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

void FileSaving::save_file_1D(double *a, const int num, const char *main_dir, const char *name, ...)
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

void FileSaving::save_file_1D(double *a, const int num, const int step, const char *main_dir, const char *name, ...)
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

void FileSaving::save_file_1D(int *a, const int num, const char *main_dir, const char *name, ...)
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

void FileSaving::save_file_1D(int *a, const int num, const int step, const char *main_dir, const char *name, ...)
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

void FileSaving::save_file_2D(double *a, int columns, int strings, const char *main_dir, const char *name, ...)
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

void FileSaving::save_file_2D_transpose(double *a, int columns, int strings, const char *main_dir, const char *name, ...)
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
