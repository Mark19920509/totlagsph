cmake_minimum_required(VERSION 2.8)

project(FDMSPH)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(PROJECT_PATH ${PROJECT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR})

set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake" ${CMAKE_MODULE_PATH})
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_FOUND)

find_package(Eigen3 REQUIRED)

IF( NOT EIGEN3_INCLUDE_DIR )
    MESSAGE( FATAL_ERROR "Please point the environment variable EIGEN3_INCLUDE_DIR to the include directory of your Eigen3 installation.")
ENDIF()

INCLUDE_DIRECTORIES ( "${EIGEN3_INCLUDE_DIR}" )

include(${PROJECT_PATH}/CMake/Common.cmake)


subdirs(src
	extern/CompactNSearch)


