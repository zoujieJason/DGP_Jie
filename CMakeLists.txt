cmake_minimum_required(VERSION 3.15)
project(DGP_Jie)

set(CMAKE_CXX_STANDARD 11)

set(DATA_PATH ${PROJECT_SOURCE_DIR}/examples/data CACHE PATH "location of shared data resources")

#Eigen
find_package(Eigen3 REQUIRED)
if ( NOT Eigen3_FOUND )
    message(STATUS "This project requires the Eigen library, and will not be compiled.")
    return()
endif()

# CGAL and its components
find_package( CGAL QUIET COMPONENTS  )
if ( NOT CGAL_FOUND )
    message(STATUS "This project requires the CGAL library, and will not be compiled.")
    return()
endif()
# Boost and its components
find_package( Boost REQUIRED )
if ( NOT Boost_FOUND )
    message(STATUS "This project requires the Boost library, and will not be compiled.")
    return()
endif()

#cinolib
include_directories(./external/include)
#add_compile_definitions(CINOLIB_USES_VTK)

#OpenMP
include_directories(/usr/local/include/)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xpreprocessor -fopenmp")

#GTest
add_subdirectory(test)

#Lib
add_subdirectory(pmp_jie)
add_subdirectory(optimizer)

#examples
add_subdirectory(examples/barycentric_mapping)
add_subdirectory(examples/lscm)
add_subdirectory(examples/abfpp)
add_subdirectory(examples/mesh_denoising)
add_subdirectory(examples/simplification)
add_subdirectory(examples/correspond)
add_subdirectory(examples/mesh_optimization_xu06)
add_subdirectory(examples/mesh_optimization_Renka16)