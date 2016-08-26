add_library(mymath mymath.cpp)

add_library(file_saving file_saving.cpp)
target_link_libraries(file_saving ${ZLIB_LIBRARIES})

add_library(pfunc pfunc.cpp)
target_link_libraries (pfunc ${PYTHON_LIBRARIES})

add_library(pyinput pyinput.cpp)
target_link_libraries(pyinput pfunc ${PYTHON_LIBRARIES})

add_library(mesh mesh.cpp)
target_link_libraries(mesh pyinput)

add_library(fdtd fdtd.cpp)
target_link_libraries(fdtd mymath pfunc pyinput mesh)

add_library(pfc pfc.cpp)
target_link_libraries(pfc mymath file_saving pfunc pyinput mesh fdtd ${ZLIB_LIBRARIES})

add_library(plasmas plasmas.cpp)
target_link_libraries(plasmas mymath pfunc pyinput mesh pfc)

add_library(particle particle.cpp)

add_library(testparticles testparticles.cpp)
target_link_libraries(testparticles pyinput mesh particle)

add_library(solver solver.cpp)
target_link_libraries(solver pyinput mesh fdtd plasmas testparticles)

add_library(output output.cpp)
add_dependencies(output pyinput file_saving mesh pfc fdtd testparticles solver)

add_executable(vlasov main.cpp)
target_link_libraries (vlasov pyinput mesh solver output)

set_property(TARGET vlasov PROPERTY RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

if (OPENMP_FOUND)
  target_compile_options(pfc PRIVATE ${OpenMP_CXX_FLAGS})
  target_compile_options(vlasov PRIVATE ${OpenMP_CXX_FLAGS})
  target_link_libraries(vlasov ${OpenMP_CXX_FLAGS})
endif (OPENMP_FOUND)
