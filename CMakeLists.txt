cmake_minimum_required(VERSION 3.16)
project(maf_cpp LANGUAGES CXX C)

# Maximum Agreement Forest solver.

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2")

# --- HiGHS 1.15.0 --------------------------------------------------------------------------------
# HiGHS 1.15.0 source is vendored in src/dependencies/HiGHS-1.15.0/.  HIPO/OpenBLAS are OFF, so the
# resulting library has NO libblastrampoline / BLAS dependency. See src/dependencies/README.md and
# THIRD_PARTY_NOTICES.md.
#
# HiGHS compiling is by far the slowest part of the build.  To avoid recompiling it on every fresh MAF
# build directory (e.g. `rm -rf build`, or a separate per-track build dir), we PREFER a PREBUILT HiGHS
# install: docker_setup.sh compiles+installs HiGHS ONCE into HIGHS_INSTALL_DIR (default ./.highs) and
# every subsequent MAF configure just LINKS it (no recompile).  Only if that prebuilt install is absent
# do we fall back to building HiGHS from the vendored source as an in-tree subproject (so a bare
# `cmake -S . -B build` still works standalone, just rebuilding HiGHS in that build tree).
set(HIGHS_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/dependencies/HiGHS-1.15.0")
set(HIGHS_INSTALL_DIR "${CMAKE_CURRENT_SOURCE_DIR}/.highs" CACHE PATH
    "Prebuilt HiGHS install prefix (built once by docker_setup.sh; linked instead of recompiled)")
if(NOT EXISTS "${HIGHS_SRC}/CMakeLists.txt")
  message(FATAL_ERROR
    "HiGHS 1.15.0 source not found at ${HIGHS_SRC} (vendored source missing from the repo).")
endif()

find_package(highs CONFIG QUIET PATHS "${HIGHS_INSTALL_DIR}" NO_DEFAULT_PATH)
if(highs_FOUND)
  message(STATUS "HiGHS: linking PREBUILT install at ${HIGHS_INSTALL_DIR} (no recompile)")
  set(MAF_HIGHS_TARGET highs::highs)
else()
  message(STATUS "HiGHS: no prebuilt install at ${HIGHS_INSTALL_DIR} -> building from vendored source "
                 "in this build tree (run docker_setup.sh to cache it in .highs and skip future recompiles)")
  # HIGHS_NO_DEFAULT_THREADS=ON builds HiGHS WITHOUT its multithreading machinery — the MAF solver is
  # single-threaded by design (matches the PACE single-thread evaluation, keeps results deterministic;
  # the solver additionally pins threads=1 at runtime).  BUILD_SHARED_LIBS=OFF + ZLIB=OFF build HiGHS as
  # a STATIC library with no zlib dependency, so the final MAF binary can be fully STATIC.
  set(BUILD_SHARED_LIBS       OFF CACHE BOOL "" FORCE)
  set(ZLIB                    OFF CACHE BOOL "" FORCE)
  set(FAST_BUILD              ON  CACHE BOOL "" FORCE)
  set(HIPO                    OFF CACHE BOOL "" FORCE)
  set(BUILD_OPENBLAS          OFF CACHE BOOL "" FORCE)
  set(BUILD_TESTING           OFF CACHE BOOL "" FORCE)
  set(BUILD_EXAMPLES          OFF CACHE BOOL "" FORCE)
  set(HIGHS_NO_DEFAULT_THREADS ON CACHE BOOL "" FORCE)
  add_subdirectory("${HIGHS_SRC}" "${CMAKE_CURRENT_BINARY_DIR}/highs_build" EXCLUDE_FROM_ALL)
  set(MAF_HIGHS_TARGET highs)
endif()

# --- core library (header-mostly modules) ------------------------------------------------------
add_library(maf_core INTERFACE)
target_include_directories(maf_core INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_link_libraries(maf_core INTERFACE ${MAF_HIGHS_TARGET})

# --- production binary --------------------------------------------------------------------------
# Track selection: -DMAF_TRACK = exact (default; prove optimality, nothing on timeout)
#                              | heuristic (5 min, best incumbent on SIGTERM)
#                              | lowerbound (10 min, integer lower bound with a,b early exit).
set(MAF_TRACK "exact" CACHE STRING "Competition track: exact, heuristic or lowerbound")
set_property(CACHE MAF_TRACK PROPERTY STRINGS exact heuristic lowerbound)

add_executable(maf "${CMAKE_CURRENT_SOURCE_DIR}/src/main.cpp")
set_target_properties(maf PROPERTIES OUTPUT_NAME "main_${MAF_TRACK}")   # binary is build/<track>/main_<track>
target_link_libraries(maf PRIVATE maf_core)
# Fully STATIC standalone binary (both tracks): no libhighs/libstdc++/libgcc runtime dependency.
target_link_options(maf PRIVATE -static -static-libstdc++ -static-libgcc)
if(MAF_TRACK STREQUAL "heuristic")
  # Always on (no build option): O(n log n) centroid 2-tree MAST (MAST-peel + branch-and-price
  # pricing, via MAF_2TREE_CENTROID) and Day linear-time reductions + cluster decomposition.
  target_compile_definitions(maf PRIVATE MAF_HEURISTIC MAF_2TREE_CENTROID)
  message(STATUS "MAF track: HEURISTIC (5 min; 2-tree MAST = centroid; linear reductions + clusters)")
elseif(MAF_TRACK STREQUAL "lowerbound")
  # Always on (no build option): centroid 2-tree MAST + Day linear-time reductions + clusters.
  target_compile_definitions(maf PRIVATE MAF_LOWERBOUND MAF_2TREE_CENTROID)
  message(STATUS "MAF track: LOWER-BOUND (10 min, approx forest; 2-tree MAST = centroid; linear reductions + clusters)")
elseif(MAF_TRACK STREQUAL "exact")
  message(STATUS "MAF track: exact (prove optimality, nothing on timeout)")
else()
  message(FATAL_ERROR "MAF_TRACK must be 'exact', 'heuristic' or 'lowerbound' (got '${MAF_TRACK}')")
endif()

# --- tests ------------------------------------------------------------------------------------
enable_testing()
file(GLOB TEST_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/tests/test_*.cpp")
foreach(tsrc ${TEST_SRCS})
  get_filename_component(tname ${tsrc} NAME_WE)
  add_executable(${tname} ${tsrc})
  target_link_libraries(${tname} PRIVATE maf_core)
  add_test(NAME ${tname} COMMAND ${tname})
endforeach()
