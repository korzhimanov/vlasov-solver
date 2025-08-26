cmake_minimum_required (VERSION 3.10)

project (Vlasov)

find_package(ZLIB REQUIRED)
find_package(Python COMPONENTS Development REQUIRED)
find_package(OpenMP)

include_directories(. ${ZLIB_INCLUDE_DIRS} ${Python_INCLUDE_DIRS})
link_directories(${Python_LIBRARY_DIRS})

add_compile_options(-Wall -g -O0 -std=c++11)

add_subdirectory(src)
