cmake_minimum_required(VERSION 3.26)
project(HW4)

set(CMAKE_CXX_STANDARD 14)

include_directories(.)

add_executable(HW4
        main.cpp
        Fock_Matrix.cpp
        Fock_Matrix.h
        Overlap_Integrals.cpp
        Overlap_Integrals.h)

include_directories(/Users/sahandadibnia/armadillo-12.6.7/include)

# Find LAPACK
find_package(LAPACK REQUIRED)

# Add Armadillo and LAPACK libraries
find_library(ARMA_LIBRARY armadillo)
target_link_libraries(HW4 ${ARMA_LIBRARY} ${LAPACK_LIBRARIES})
