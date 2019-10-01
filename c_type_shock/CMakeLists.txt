# -----------------------------------------------------------------
# $Revision: 4483 $
# $Date: 2015-04-22 16:30:17 -0700 (Wed, 22 Apr 2015) $
# -----------------------------------------------------------------

cmake_minimum_required(VERSION 2.8)

# Specify project name
PROJECT(c_type_shock CXX)

MARK_AS_ADVANCED(EXECUTABLE_OUTPUT_PATH LIBRARY_OUTPUT_PATH)

# Set names of examples specific dependency source files
SET(examples_dependencies chemistry.cpp chemical_reaction_data.cpp coll_rates_v2.cpp coll_rates_h2_v2.cpp coll_rates_h2o_v2.cpp coll_rates_ions.cpp coll_rates_co.cpp differential_equations.cpp dust_model_v2.cpp photoelectric_emission.cpp elastic_scattering.cpp lvg_method_functions.cpp radiation_field.cpp spectroscopy_v2.cpp integration.cpp interpolation.cpp special_functions.cpp linear_algebra.cpp reaction_rate_analyzer.cpp)

# Specify path to SUNDIALS header files
SET(SUNDIALS_INC_DIR
  /usr_homes/alexn/Sundials/cvode-2.8.2-install/include
  CACHE STRING
  "Location of SUNDIALS header files")

# Add path to SUNDIALS header files
INCLUDE_DIRECTORIES(${SUNDIALS_INC_DIR})

# Set search path for SUNDIALS libraries 
SET(SUNDIALS_LIB_DIR /usr_homes/alexn/Sundials/cvode-2.8.2-install/lib)

# Find the SUNDIALS solver's library
FIND_LIBRARY(SUNDIALS_SOLVER_LIB
  sundials_cvode ${SUNDIALS_LIB_DIR}
  DOC "CVODE library")

# Find the NVECTOR library
FIND_LIBRARY(SUNDIALS_NVEC_LIB
  sundials_nvecserial ${SUNDIALS_LIB_DIR}
  DOC "NVECTOR library")

# Set additional libraries
SET(SUNDIALS_EXTRA_LIB  -lm /usr/lib64/librt.so CACHE STRING "Additional libraries")

# Set other libraries libraries
SET(LAPACK_LIBRARIES  CACHE STRING "Lapack libraries")
SET(SUPERLUMT_LIBRARIES  CACHE STRING "SuperLUMT libraries")
SET(KLU_LIBRARIES  CACHE STRING "KLU libraries")

# List of all libraries
SET(SUNDIALS_LIBS ${SUNDIALS_SOLVER_LIB} ${SUNDIALS_NVEC_LIB} ${SUNDIALS_EXTRA_LIB})
IF(LAPACK_LIBRARIES)
  LIST(APPEND SUNDIALS_LIBS ${LAPACK_LIBRARIES})
ENDIF(LAPACK_LIBRARIES)

IF(SUPERLUMT_LIBRARIES)
  LIST(APPEND SUNDIALS_LIBS ${SUPERLUMT_LIBRARIES})
ENDIF(SUPERLUMT_LIBRARIES)

IF(KLU_LIBRARIES)
  LIST(APPEND SUNDIALS_LIBS ${KLU_LIBRARIES})
ENDIF(KLU_LIBRARIES)

find_package(OpenMP)
if(OPENMP_FOUND)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

ADD_EXECUTABLE(c_type_shock c_type_shock.cpp ${examples_dependencies})
TARGET_LINK_LIBRARIES(c_type_shock ${SUNDIALS_LIBS})
