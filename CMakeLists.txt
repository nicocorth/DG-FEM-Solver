PROJECT(GMSHAPI CXX)
CMAKE_MINIMUM_REQUIRED(VERSION 3.1)

# put all generated files in a single folder (required in windows for exe linked to dlls)
SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "")
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin CACHE PATH "")
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

# build type is "" by default in Linux
IF(NOT CMAKE_BUILD_TYPE)
    SET( CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
ENDIF()

# enable C++11
SET(CMAKE_CXX_STANDARD 11) # newer way to set C++11 (requires cmake>3.1)
SET(CMAKE_CXX_STANDARD_REQUIRED ON)

# compiler options
IF(UNIX AND CMAKE_CXX_COMPILER_ID STREQUAL GNU)
    # allows gcc5/ubuntu to link against the c++ interface
    # (this does not work with mingw)
    ADD_DEFINITIONS(-D_GLIBCXX_USE_CXX11_ABI=1) # libgmsh is linked to libgcc4 (pbl with std::string) 
ELSEIF(CMAKE_CXX_COMPILER_ID MATCHES MSVC)
    ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS -D_SCL_SECURE_NO_DEPRECATE -D_CRT_NONSTDC_NO_WARNINGS)
    ADD_DEFINITIONS(-D_USE_MATH_DEFINES) # for M_PI
ENDIF()

IF(APPLE)
    # on macOS, do not give priority to frameworks/apps
    SET(CMAKE_FIND_APPBUNDLE LAST)
    SET(CMAKE_FIND_FRAMEWORK LAST)
ENDIF()

message(STATUS "Checking OpenMP")
find_package(OpenMP)
IF(OPENMP_FOUND)
    message("OpenMP found.")
    # add flags for OpenMP
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${OpenMP_SHARED_LINKER_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
ELSE()
    message("OpenMP missing")
ENDIF()

# find gmsh-sdk
# gmsh.h
INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/gmsh-4.4.1-Linux64-sdk/include")
# headers
INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/Headers")
# eigen include files
INCLUDE_DIRECTORIES("${CMAKE_CURRENT_SOURCE_DIR}/eigen")

FIND_LIBRARY(GMSH_LIBRARIES "libgmsh.so.4.4" PATH "${CMAKE_CURRENT_SOURCE_DIR}/gmsh-4.4.1-Linux64-sdk/lib")
MESSAGE(STATUS "GMSH_LIBRARIES=" ${GMSH_LIBRARIES})
IF(NOT GMSH_LIBRARIES)
    MESSAGE(FATAL_ERROR "gmsh library 4.4.1 not found!")
ENDIF()

# find the gmsh executable
FIND_PROGRAM(GMSH_EXECUTABLE gmsh PATH "${CMAKE_CURRENT_SOURCE_DIR}/gmsh-4.4.1-Linux64-sdk/bin")
MESSAGE(STATUS "GMSH_EXECUTABLE=" ${GMSH_EXECUTABLE})
IF(NOT GMSH_EXECUTABLE)
    MESSAGE(FATAL_ERROR "gmsh executable not found!")
ENDIF()


FILE(GLOB SRCS "main.cpp" "Mesh/*.cpp" "Solver/*.cpp")
SET(PROG_NAME main)
ADD_EXECUTABLE(${PROG_NAME} ${SRCS})
TARGET_LINK_LIBRARIES(${PROG_NAME} ${GMSH_LIBRARIES})
