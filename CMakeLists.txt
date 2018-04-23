# CMake compatibility issues: don't modify this, please!
cmake_minimum_required(VERSION 2.6.3)

project(brownie)

# project version
set(${PROJECT_NAME}_MAJOR_VERSION 0)
set(${PROJECT_NAME}_MINOR_VERSION 1)
set(${PROJECT_NAME}_PATCH_LEVEL 0)

# set the default configuration to Release
if(NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release" CACHE STRING
            "Choose the type of build, options are: Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE)

# set the module path
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

# set some definitions
if(MAXKMERLENGTH)
        add_definitions("-DMAXKMERLENGTH=${MAXKMERLENGTH}")
else(MAXKMERLENGTH)
        add_definitions("-DMAXKMERLENGTH=63")
endif(MAXKMERLENGTH)
add_definitions("-DBROWNIE_MAJOR_VERSION=${${PROJECT_NAME}_MAJOR_VERSION}")
add_definitions("-DBROWNIE_MINOR_VERSION=${${PROJECT_NAME}_MINOR_VERSION}")
add_definitions("-DBROWNIE_PATCH_LEVEL=${${PROJECT_NAME}_PATCH_LEVEL}")

# set windows specific flags
if (MSVC)
        add_definitions("-D_SCL_SECURE_NO_WARNINGS")
        add_definitions("-D_CRT_SECURE_NO_WARNINGS")
endif (MSVC)

# add_definitions(-DDEBUG)        # FIXME: remove !!

# set g++ specific flags
if (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)
        set(CMAKE_CXX_FLAGS "-Wno-deprecated -std=c++11 -fopenmp")
        set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3 -Wall -pedantic -Wno-long-long")
        set(CMAKE_CXX_FLAGS_RELEASE "-O3")
        set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g3")
        set(CMAKE_C_FLAGS_DEBUG "-O0 -g3 -Wall -pedantic -Wno-long-long")
endif (CMAKE_COMPILER_IS_GNUCXX OR CMAKE_COMPILER_IS_GNUCC)

# uncomment the portion below to disable assertions
if (CMAKE_BUILD_TYPE STREQUAL Release)
        add_definitions(-DNDEBUG)
else (CMAKE_BUILD_TYPE STREQUAL Release)
        add_definitions(-DDEBUG)
endif (CMAKE_BUILD_TYPE STREQUAL Release)

# check if zlib is present
find_package(ZLIB)
if (ZLIB_FOUND)
    add_definitions(-DHAVE_ZLIB)
    include_directories(${ZLIB_INCLUDE_DIR})
endif(ZLIB_FOUND)

set(CMAKE_VERBOSE_MAKEFILE ON)

# set include path for Google's sparse hash table
find_package(SparseHash)
if (SPARSEHASH_FOUND)
    include_directories(${SPARSEHASH_INCLUDE_DIR})
else (SPARSEHASH_FOUND)
    message(FATAL_ERROR "\nFATAL ERROR: The required Google SparseHash package"
            " could not be found on this system.  Please refer to the Velvet "
            "manual for the Google Sparsehash installation instructions.  If "
            "you installed Google Sparsehash in a non-standard location "
            "(e.g. somewhere in your homedir), you can point cmake to the "
            "installation location as follows: \ncmake "
            "-DSPARSEHASH_INCLUDE_DIR=<path-to-sparsehash>/include .")
endif(SPARSEHASH_FOUND)

add_subdirectory(src)
add_subdirectory(unittest)
