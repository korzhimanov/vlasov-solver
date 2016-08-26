cmake_minimum_required (VERSION 2.8)

project (Vlasov)

find_package(ZLIB REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(OpenMP)

include_directories(./include ${ZLIB_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS})

add_compile_options(-Wall -g -O0)

add_subdirectory(src)
