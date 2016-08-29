/*
 * Use of this source code is governed by the MIT License (MIT).
 * See the LICENSE file for details.
 *
 * Copyright (c) 2016 Artem Korzhimanov <korzhimanov.artem@gmail.com>
 */

/**
 * \file errors.h
 * \brief The header file which defines enum with errors codes
 * \author Artem Korzhimanov
 * \copyright The MIT License (MIT)
 */

#ifndef INCLUDE_ERRORS_H_
#define INCLUDE_ERRORS_H_

typedef enum {
  // pyinput
  PYTHON_INITIALIZATION_FAILED,
  CANNOT_OPEN_INIT_FILE,
  VAR_NOT_FOUND,
  VAR_MUST_INTEGER,
  VAR_MUST_FLOAT,
  VAR_MUST_POSITIVE,
  VAR_MUST_NOTNEGATIVE,
  VAR_MUST_BOOLEAN,
  VAR_MUST_STRING,
  VAR_MUST_FUNCTION,

  // Solver
  THETA_OUT_OF_RANGE,

  // Output
  WRONG_OUTPUT_FORMAT
} ErrorCodes;

#endif  // INCLUDE_ERRORS_H_
