cmake_minimum_required (VERSION 2.6)
project (spectral_scalar_toy)

set(LIBRARY ScalarToyLib)

set(LIBRARY_SOURCES
    scalarFunction)

set(SYSTEM_LINKS
    -lboost_iostreams
    -lboost_system
    -lboost_filesystem)

add_library(${LIBRARY} ${LIBRARY_SOURCES})


target_link_libraries(${LIBRARY} ${SYSTEM_LINKS})

set(LIBS_TO_LINK ScalarToyLib)

add_executable(run_wave run_wave.cpp)

target_link_libraries(run_wave ${LIBS_TO_LINK})

set(CMAKE_BUILD_TYPE Debug)
