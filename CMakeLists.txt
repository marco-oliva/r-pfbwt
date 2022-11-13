cmake_minimum_required(VERSION 3.15)

project(r-pfbwt)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

################################################################################
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose Release or Debug" FORCE)
endif()
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

################################################################################
# Version number
SET(VERSION_MAJOR "0")
SET(VERSION_MINOR "1")
SET(VERSION_PATCH "0")
SET(VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
set(PROJECT_VERSION ${VERSION})

execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE)

execute_process(
        COMMAND git log -1 --format=%h
        WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE)

message(STATUS "Git current branch: ${GIT_BRANCH}")
message(STATUS "Git commit hash: ${GIT_COMMIT_HASH}")

message(STATUS "Generating version.hpp")
configure_file(
        ${PROJECT_SOURCE_DIR}/include/version.hpp.in
        ${PROJECT_BINARY_DIR}/generated/version.hpp)

include_directories(${PROJECT_BINARY_DIR}/generated)

################################################################################
# Compiler and linker flags

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    set(BUILD_TYPE_FLAGS "-O3")
elseif(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(BUILD_TYPE_FLAGS "-g")
elseif(CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(BUILD_TYPE_FLAGS "-O3 -g -DNDEBUG")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${BUILD_TYPE_FLAGS}")


################################################################################
# Thirdparty

include(FetchContent)

# gsacak
FetchContent_Declare(
        gsacak
        GIT_REPOSITORY https://github.com/felipelouza/gsa-is.git
)
FetchContent_GetProperties(gsacak)
if(NOT gsacak_POPULATED)
    FetchContent_Populate(gsacak)
    add_library(gsacak OBJECT ${gsacak_SOURCE_DIR}/gsacak.c ${gsacak_SOURCE_DIR}/gsacak.h)
    target_include_directories(gsacak PUBLIC "${gsacak_SOURCE_DIR}")

    add_library(gsacak64 OBJECT ${gsacak_SOURCE_DIR}/gsacak.c ${gsacak_SOURCE_DIR}/gsacak.h)
    target_include_directories(gsacak64 PUBLIC "${gsacak_SOURCE_DIR}")
    target_compile_options(gsacak64 PUBLIC -DM64)
endif()

# sdsl
FetchContent_Declare(
        sdsl
        GIT_REPOSITORY https://github.com/simongog/sdsl-lite
)
FetchContent_GetProperties(sdsl)
if(NOT sdsl_POPULATED)
    FetchContent_Populate(sdsl)
    set(GENERATE_DOC OFF CACHE BOOL "Do not generate doxygen for sdsl-lite")
    add_subdirectory(${sdsl_SOURCE_DIR} ${sdsl_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()

# divsuffsort
FetchContent_Declare(
        divsufsort
        GIT_REPOSITORY https://github.com/simongog/libdivsufsort.git
        GIT_TAG        2.0.1
)
FetchContent_GetProperties(divsufsort)
if(NOT divsufsort_POPULATED)
    FetchContent_Populate(divsufsort)
    set(BUILD_SHARED_LIBS OFF CACHE BOOL "Do not build a shared library for libdivsufsort")
    set(BUILD_EXAMPLES OFF CACHE BOOL "Do not build libdivsufsort example")
    set(BUILD_DIVSUFSORT64 ON CACHE BOOL "Build libdivsufsort in 64-bits mode")
    add_subdirectory(${divsufsort_SOURCE_DIR} ${divsufsort_BINARY_DIR} EXCLUDE_FROM_ALL)
    target_include_directories(divsufsort PUBLIC "${divsufsort_BINARY_DIR}/include")
    target_include_directories(divsufsort64 PUBLIC "${divsufsort_BINARY_DIR}/include")
endif()

# malloc count
FetchContent_Declare(
        malloc_count
        GIT_REPOSITORY https://github.com/bingmann/malloc_count.git
        GIT_TAG        v0.7.1
)
FetchContent_GetProperties(malloc_count)
if(NOT malloc_count_POPULATED)
    FetchContent_Populate(malloc_count)
    add_library(malloc_count OBJECT ${malloc_count_SOURCE_DIR}/malloc_count.c ${malloc_count_SOURCE_DIR}/malloc_count.h)
    target_include_directories(malloc_count PUBLIC "${malloc_count_SOURCE_DIR}")
    target_link_libraries(malloc_count dl)
endif()

# Catch
FetchContent_Declare(
        Catch2
        GIT_REPOSITORY https://github.com/catchorg/Catch2.git
        GIT_TAG        v3.0.1
)
# Making available only if building tests
FetchContent_GetProperties(Catch2)
if (NOT Catch2_POPULATED AND (BUILD_TESTS OR ${CMAKE_BUILD_TYPE} MATCHES "Debug"))
    FetchContent_MakeAvailable(Catch2)
endif()

# SPDLOG
FetchContent_Declare(
        spdlog
        GIT_REPOSITORY https://github.com/gabime/spdlog.git
        GIT_TAG        v1.10.0
)
FetchContent_GetProperties(spdlog)
if(NOT spdlog_POPULATED)
    FetchContent_Populate(spdlog)
endif()
include_directories(${spdlog_SOURCE_DIR}/include)

# CLI11
FetchContent_Declare(
        cli11
        GIT_REPOSITORY "https://github.com/CLIUtils/CLI11.git"
        GIT_TAG "v1.9.0"
)
FetchContent_GetProperties(cli11)
if(NOT cli11_POPULATED)
    FetchContent_Populate(cli11)
endif()
include_directories(${cli11_SOURCE_DIR}/include)

# PFP-DS
FetchContent_Declare(
        pfp_ds
        GIT_REPOSITORY git@github.com:marco-oliva/PFP-DS.git
        GIT_TAG        easy_hard_algo
)
FetchContent_GetProperties(pfp_ds)
if(NOT pfp_ds_POPULATED)
    FetchContent_Populate(pfp_ds)
endif()
include_directories(${pfp_ds_SOURCE_DIR}/include)


# mio mmap
FetchContent_Declare(
        mio
        GIT_REPOSITORY https://github.com/mandreyel/mio.git
)
FetchContent_GetProperties(mio)
if(NOT mio_POPULATED)
    FetchContent_Populate(mio)
endif()
include_directories(${mio_SOURCE_DIR}/single_include)


# rlestring
FetchContent_Declare(
        rlestring
        GIT_REPOSITORY https://github.com/marco-oliva/rle_string.git
        GIT_TAG        main
)
FetchContent_GetProperties(rlestring)
if(NOT rlestring_POPULATED)
    FetchContent_Populate(rlestring)
    add_subdirectory(${rlestring_SOURCE_DIR} ${rlestring_BINARY_DIR} EXCLUDE_FROM_ALL)
endif()
include_directories(${rlestring_SOURCE_DIR}/include)

# OpenMP
find_package(OpenMP)
if(OPENMP_FOUND)
    include_directories(${OpenMP_CXX_INCLUDE_DIRS})
    link_libraries(${OpenMP_CXX_LIBRARIES})
    if(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    elseif(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        # assumes clang build
        # we can't reliably detect when we're using clang, so for the time being we assume
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -Xpreprocessor -fopenmp")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -Xpreprocessor -fopenmp")
        set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif()
endif()

FetchContent_Declare(
        kseq
        GIT_REPOSITORY https://github.com/lh3/seqtk.git
        GIT_TAG        v1.3
)
if(NOT kseq_POPULATED)
    FetchContent_Populate(kseq)
endif()
include_directories(${kseq_SOURCE_DIR})

################################################################################
# Lib

include_directories(${PROJECT_SOURCE_DIR}/include)
include_directories(${PROJECT_BINARY_DIR}/generated)

file(GLOB SOURCES ${CMAKE_SOURCE_DIR}/src/*.cpp)
if (ENABLE_MALLOC_COUNT)
    message(STATUS "Using malloc_count.c")
    file(GLOB MALLOC_COUNT_SOURCES ${mcount_SOURCE_DIR}/*.c)
else()
    set(MALLOC_COUNT_SOURCES "")
endif()

# add_library(rpfbwt ${SOURCES} ${MALLOC_COUNT_SOURCES} ${MURMUR_SOURCES})

################################################################################
# Tests

if (BUILD_TESTS OR ${CMAKE_BUILD_TYPE} MATCHES "Debug")
    message(STATUS "Building unit tests")

    add_executable(unit_tests_rpfbwt tests/unit_tests.cpp)
    target_link_libraries(unit_tests_rpfbwt Catch2::Catch2 sdsl gsacak divsufsort divsufsort64 rlestring)

    add_executable(unit_tests_rpfbwt_64 tests/unit_tests.cpp)
    target_link_libraries(unit_tests_rpfbwt_64 Catch2::Catch2 sdsl gsacak64 divsufsort divsufsort64 rlestring)

    enable_testing()
    add_test(NAME unit_tests_rpfbwt COMMAND unit_tests_rpfbwt)
endif()

################################################################################
# Targets

add_executable(rpfbwt32 rpfbwt.cpp)
target_link_libraries(rpfbwt32 dl gsacak sdsl divsufsort divsufsort64 rlestring)

add_executable(rpfbwt64 rpfbwt.cpp)
target_link_libraries(rpfbwt64 dl gsacak64 sdsl divsufsort divsufsort64 rlestring)

################################################################################
# Install

install(TARGETS rpfbwt32 rpfbwt64)

################################################################################
# CPack

set(CPACK_PACKAGE_VERSION_MAJOR ${VERSION_MAJOR})
set(CPACK_PACKAGE_VERSION_MINOR ${VERSION_MINOR})
set(CPACK_PACKAGE_VERSION_PATCH ${VERSION_PATCH})
set(CPACK_PACKAGE_VERSION "${VERSION}")

include(InstallRequiredSystemLibraries)
set(CPACK_GENERATOR "STGZ;TGZ;DEB")
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_PACKAGE_VENDOR "University of Florida")
set(CPACK_PACKAGE_CONTACT "marco.oliva@ufl.edu")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Big-BWT")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE")
set(CPACK_RESOURCE_FILE_README "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_PACKAGE_NAME "${CMAKE_PROJECT_NAME}")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-sources")

set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Marco Oliva")
set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)
set(CPACK_COMPONENTS_GROUPING ALL_COMPONENTS_IN_ONE)
set (CPACK_DEBIAN_PACKAGE_SHLIBDEPS ON)
set(CPACK_DEB_COMPONENT_INSTALL YES)
include(CPack)





