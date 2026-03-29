cmake_minimum_required(VERSION 3.16)
project(linear_smoothing_sfem_cpp LANGUAGES CXX)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(PkgConfig QUIET)
include(CheckIncludeFileCXX)
get_target_property(EIGEN3_TARGET_INCLUDE_DIRS Eigen3::Eigen INTERFACE_INCLUDE_DIRECTORIES)
option(ENABLE_OPENMP "Enable OpenMP parallel assembly when available" ON)
option(ENABLE_UMFPACK "Enable SuiteSparse/UMFPACK when available" ON)

if(ENABLE_UMFPACK AND PkgConfig_FOUND)
    pkg_check_modules(UMFPACK_PC QUIET UMFPACK)
    if(NOT UMFPACK_PC_FOUND)
        pkg_check_modules(UMFPACK_PC QUIET umfpack)
    endif()
endif()

if(ENABLE_UMFPACK)
    if(UMFPACK_PC_FOUND)
        set(UMFPACK_FOUND TRUE)
        set(UMFPACK_INCLUDE_DIRS ${UMFPACK_PC_INCLUDE_DIRS})
        set(UMFPACK_LINK_LIBRARIES ${UMFPACK_PC_LINK_LIBRARIES})
        set(UMFPACK_DETECTION_METHOD "pkg-config")
    else()
        find_path(UMFPACK_INCLUDE_DIR
            NAMES umfpack.h
            HINTS
                /opt/homebrew/include
                /opt/homebrew/include/suitesparse
                /usr/local/include
                /usr/local/include/suitesparse
                /usr/include
                /usr/include/suitesparse
            PATH_SUFFIXES suitesparse
        )
        find_library(UMFPACK_LIBRARY
            NAMES umfpack
            HINTS
                /opt/homebrew/lib
                /usr/local/lib
                /usr/lib
        )

        if(UMFPACK_INCLUDE_DIR AND UMFPACK_LIBRARY)
            set(UMFPACK_FOUND TRUE)
            set(UMFPACK_INCLUDE_DIRS ${UMFPACK_INCLUDE_DIR})
            set(UMFPACK_LINK_LIBRARIES ${UMFPACK_LIBRARY})
            set(UMFPACK_DETECTION_METHOD "manual search")
        endif()
    endif()

    set(CMAKE_REQUIRED_INCLUDES ${EIGEN3_TARGET_INCLUDE_DIRS})
    if(UMFPACK_FOUND)
        list(APPEND CMAKE_REQUIRED_INCLUDES ${UMFPACK_INCLUDE_DIRS})
    endif()
    check_include_file_cxx("Eigen/UmfPackSupport" EIGEN_HAS_UMFPACK_SUPPORT_NEW)
    if(NOT EIGEN_HAS_UMFPACK_SUPPORT_NEW)
        check_include_file_cxx("unsupported/Eigen/UmfPackSupport" EIGEN_HAS_UMFPACK_SUPPORT_OLD)
    endif()
    unset(CMAKE_REQUIRED_INCLUDES)

    if(EIGEN_HAS_UMFPACK_SUPPORT_NEW OR EIGEN_HAS_UMFPACK_SUPPORT_OLD)
        set(EIGEN_HAS_UMFPACK_SUPPORT TRUE)
    else()
        set(EIGEN_HAS_UMFPACK_SUPPORT FALSE)
    endif()
endif()

if(ENABLE_UMFPACK AND UMFPACK_FOUND AND EIGEN_HAS_UMFPACK_SUPPORT)
    add_library(SuiteSparse::UMFPACK INTERFACE IMPORTED)
    set_target_properties(SuiteSparse::UMFPACK PROPERTIES
        INTERFACE_INCLUDE_DIRECTORIES "${UMFPACK_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${UMFPACK_LINK_LIBRARIES}"
        INTERFACE_COMPILE_DEFINITIONS "USE_EIGEN_UMFPACK"
    )
    message(STATUS "UMFPACK enabled via ${UMFPACK_DETECTION_METHOD}")
elseif(ENABLE_UMFPACK)
    message(STATUS "UMFPACK requested but unavailable in this Eigen/toolchain setup; building without UMFPACK")
endif()

if(ENABLE_OPENMP)
    find_package(OpenMP QUIET)
endif()

if(ENABLE_OPENMP AND NOT OpenMP_CXX_FOUND AND APPLE)
    find_path(OPENMP_FALLBACK_INCLUDE_DIR
        NAMES omp.h
        HINTS
            /opt/homebrew/opt/libomp/include
            /usr/local/opt/libomp/include
    )
    find_library(OPENMP_FALLBACK_LIBRARY
        NAMES omp gomp iomp5
        HINTS
            /opt/homebrew/opt/libomp/lib
            /usr/local/opt/libomp/lib
    )

    if(OPENMP_FALLBACK_INCLUDE_DIR AND OPENMP_FALLBACK_LIBRARY)
        add_library(OpenMP::OpenMP_CXX INTERFACE IMPORTED)
        if(CMAKE_CXX_COMPILER_ID MATCHES "AppleClang|Clang")
            set(_openmp_compile_options "-Xpreprocessor;-fopenmp")
        elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            set(_openmp_compile_options "-fopenmp")
        else()
            set(_openmp_compile_options "")
        endif()

        set_target_properties(OpenMP::OpenMP_CXX PROPERTIES
            INTERFACE_COMPILE_OPTIONS "${_openmp_compile_options}"
            INTERFACE_INCLUDE_DIRECTORIES "${OPENMP_FALLBACK_INCLUDE_DIR}"
            INTERFACE_LINK_LIBRARIES "${OPENMP_FALLBACK_LIBRARY}"
        )
        set(OpenMP_CXX_FOUND TRUE)
        message(STATUS "Using OpenMP fallback library: ${OPENMP_FALLBACK_LIBRARY}")
    endif()
endif()

add_executable(csfem
    main.cpp
    src/structural_problem_library.cpp
    src/structural_solver.cpp
    src/sfem_basis.cpp
    src/sfem_edge_assembler.cpp
    src/sfem_framework.cpp
    src/sfem_io.cpp
    src/sfem_material.cpp
    src/sfem_mesh.cpp
    src/sfem_solver_utils.cpp
    src/vtk_writer.cpp
)

add_executable(esfem
    main.cpp
    src/structural_problem_library.cpp
    src/structural_solver.cpp
    src/sfem_basis.cpp
    src/sfem_edge_assembler.cpp
    src/sfem_framework.cpp
    src/sfem_io.cpp
    src/sfem_material.cpp
    src/sfem_mesh.cpp
    src/sfem_solver_utils.cpp
    src/vtk_writer.cpp
)
target_include_directories(csfem PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(csfem PRIVATE
    Eigen3::Eigen
)

target_compile_definitions(csfem PRIVATE
    FIXED_METHOD_CSFEM
)

target_include_directories(esfem PRIVATE
    ${PROJECT_SOURCE_DIR}/include
)

target_link_libraries(esfem PRIVATE
    Eigen3::Eigen
)

target_compile_definitions(esfem PRIVATE
    FIXED_METHOD_ESFEM
)

if(UMFPACK_FOUND AND EIGEN_HAS_UMFPACK_SUPPORT)
    target_link_libraries(csfem PRIVATE SuiteSparse::UMFPACK)
    target_link_libraries(esfem PRIVATE SuiteSparse::UMFPACK)
endif()

if(OpenMP_CXX_FOUND)
    target_link_libraries(csfem PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(esfem PRIVATE OpenMP::OpenMP_CXX)
    message(STATUS "OpenMP enabled")
elseif(ENABLE_OPENMP)
    message(STATUS "OpenMP requested but not found; building without OpenMP")
endif()

set_target_properties(csfem PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin
)

set_target_properties(esfem PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin
)
