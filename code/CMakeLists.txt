CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(Aspect REQUIRED 
  HINTS ${ASPECT_DIR} ../../eclipse  $ENV{ASPECT_DIR}
  )

DEAL_II_INITIALIZE_CACHED_VARIABLES()
PROJECT(aegean)
ADD_LIBRARY(aegean SHARED table_vp.cc multicomponent_vp.cc
adiabat_ascii.cc ascii_T_regions.cc ascii_data_function.cc
two_merged_chunks.cc 
point_values_eigen_vectors.cc)
ASPECT_SETUP_PLUGIN(aegean)


