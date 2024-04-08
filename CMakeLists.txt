cmake_minimum_required(VERSION 3.19)

# disallow building in same directory as source.
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_BINARY_DIR})
  message(WARNING "ERROR: CMake build directory cannot be the same as root source directory. HINT: Can build with `make`.")
endif()

set(CMAKE_DISABLE_SOURCE_CHANGES ON)
set(CMAKE_DISABLE_IN_SOURCE_BUILD ON)

project(bito)
find_package(Python COMPONENTS Interpreter Development)
set(BITO_VERSION "0.1")

set(CMAKE_SKIP_RPATH ON)
execute_process(COMMAND ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/get_rpath.py
  OUTPUT_VARIABLE BITO_RPATH
  COMMAND_ERROR_IS_FATAL ANY
  OUTPUT_STRIP_TRAILING_WHITESPACE)

include(ExternalProject)
ExternalProject_Add(beagle-lib
  GIT_REPOSITORY https://github.com/beagle-dev/beagle-lib.git
  GIT_TAG origin/hmc-clock
  GIT_SHALLOW true
  GIT_PROGRESS true
  UPDATE_DISCONNECTED true
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/beagle-lib/install -DBUILD_CUDA=OFF -DBUILD_OPENCL=OFF -DBUILD_JNI=OFF -DCMAKE_CXX_FLAGS="-DHAVE_CPUID_H" -DCMAKE_SHARED_LINKER_FLAGS="${BITO_RPATH}"
  PREFIX beagle-lib
  INSTALL_DIR beagle-lib/install
)

ExternalProject_Add(pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG origin/v2.7
  GIT_SHALLOW true
  GIT_PROGRESS true
  UPDATE_DISCONNECTED true
  CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${PROJECT_BINARY_DIR}/pybind11/install -DPYBIND11_TEST=OFF
  PREFIX pybind11
  INSTALL_DIR pybind11/install
)

option(WERROR "Treat warnings as errors" ON)

function(bito_compile_opts PRODUCT)
  target_compile_features(${PRODUCT} PUBLIC cxx_std_17)
  target_compile_options(${PRODUCT} PUBLIC -Wno-unknown-warning -Wno-unknown-warning-option -Wall)

  if(${WERROR})
    target_compile_options(${PRODUCT} PUBLIC -Werror)
  else()
    target_compile_options(${PRODUCT} PUBLIC -Wold-style-cast)
  endif()

  target_include_directories(${PRODUCT} PUBLIC
    ${PROJECT_BINARY_DIR}/beagle-lib/install/include/libhmsbeagle-1
    lib/eigen
  )
endfunction()

function(bito_link_opts PRODUCT)
  target_link_options(${PRODUCT} PUBLIC -pthread ${BITO_RPATH})
  target_link_libraries(${PRODUCT} PUBLIC bito-core)
  add_dependencies(${PRODUCT} bito-core)
endfunction()

function(bito_executable PRODUCT)
  add_executable(${PRODUCT} ${ARGN})
  bito_compile_opts(${PRODUCT})
  bito_link_opts(${PRODUCT})
endfunction()

function(bito_extra PRODUCT)
  bito_executable(${PRODUCT} ${ARGN})
  target_include_directories(${PRODUCT} PUBLIC ../src)
endfunction()

# ##################
# libbito-core.so #
# ##################
add_library(bito-core SHARED
  src/alignment.cpp
  src/bitset.cpp
  src/block_model.cpp
  src/block_specification.cpp
  src/choice_map.cpp
  src/clock_model.cpp
  src/combinatorics.cpp
  src/csv.cpp
  src/driver.cpp
  src/engine.cpp
  src/fat_beagle.cpp
  src/gp_dag.cpp
  src/gp_engine.cpp
  src/gp_instance.cpp
  src/gp_operation.cpp
  src/graft_dag.cpp
  src/mersenne_twister.cpp
  src/node.cpp
  src/numerical_utils.cpp
  src/nni_engine.cpp
  src/nni_operation.cpp
  src/parser.cpp
  src/phylo_flags.cpp
  src/phylo_model.cpp
  src/pv_handler.cpp
  src/psp_indexer.cpp
  src/quartet_hybrid_request.cpp
  src/reindexer.cpp
  src/sankoff_handler.cpp
  src/rooted_gradient_transforms.cpp
  src/rooted_sbn_instance.cpp
  src/rooted_tree.cpp
  src/rooted_tree_collection.cpp
  src/sbn_maps.cpp
  src/sbn_probability.cpp
  src/sbn_support.cpp
  src/scanner.cpp
  src/site_model.cpp
  src/site_pattern.cpp
  src/stick_breaking_transform.cpp
  src/subsplit_dag.cpp
  src/substitution_model.cpp
  src/taxon_name_munging.cpp
  src/tidy_subsplit_dag.cpp
  src/topology_sampler.cpp
  src/tp_engine.cpp
  src/tree.cpp
  src/tree_collection.cpp
  src/unrooted_sbn_instance.cpp
  src/unrooted_tree.cpp
  src/unrooted_tree_collection.cpp
  src/zlib_stream.cpp
)
bito_compile_opts(bito-core)
target_link_options(bito-core PUBLIC -pthread ${BITO_RPATH})
target_link_libraries(bito-core PUBLIC hmsbeagle hmsbeagle-cpu z)
target_link_directories(bito-core PUBLIC
  ${PROJECT_BINARY_DIR}/beagle-lib/install/lib
)
add_dependencies(bito-core beagle-lib)

# #############
# libbito.so #
# #############
add_library(bito SHARED
  src/pybito.cpp
)
bito_compile_opts(bito)
target_include_directories(bito PUBLIC
  ${Python_INCLUDE_DIRS}
  ${PROJECT_BINARY_DIR}/pybind11/install/include
)
bito_link_opts(bito)

if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
  target_link_options(bito PUBLIC -undefined dynamic_lookup)
endif()

add_dependencies(bito pybind11)

# ##############
# executables #
# ##############
bito_executable(doctest
  src/doctest.cpp)

bito_executable(gp_doctest
  src/gp_doctest.cpp)

# ##################
# optional extras #
# ##################
add_subdirectory(extras)

file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/py)
set(BITO_BINARY_DIR ${PROJECT_BINARY_DIR})
configure_file(setup.py.in
  ${PROJECT_BINARY_DIR}/setup.py @ONLY
)
configure_file(__init__.py.in
  ${PROJECT_BINARY_DIR}/py/__init__.py @ONLY
)

add_custom_target(pip ALL
  COMMAND ln -sf ${PROJECT_BINARY_DIR}/beagle-lib/install/lib/* ${PROJECT_BINARY_DIR}/py
  COMMAND ln -sf ${PROJECT_BINARY_DIR}/$<TARGET_FILE_NAME:bito-core> ${PROJECT_BINARY_DIR}/py
  COMMAND ln -sf ${PROJECT_BINARY_DIR}/$<TARGET_FILE_NAME:bito> ${PROJECT_BINARY_DIR}/py/__init__.so
  DEPENDS bito
)
