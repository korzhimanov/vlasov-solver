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
} // namespace filesaving
