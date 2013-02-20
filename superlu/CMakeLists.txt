#CMakeLists configuration file for SuperLU. Put this file in the SRC directory of SUPERLU's source tree. Configure and generate using your favourite build system.
#Author: Adam Harrison - included in LibMIA for conveinence
cmake_minimum_required (VERSION 2.8)
project(SuperLU)

#ask for location of blas lib
SET(blas_LIB  CACHE FILEPATH /path/to/blas_lib)

#get all .h and .c files
file(GLOB SuperLU_SRC
    "*.h"
    "*.c"
)
#create a static target and link to blas
add_library(SuperLU STATIC ${SuperLU_SRC})
target_link_libraries(SuperLU ${blas_LIB})


