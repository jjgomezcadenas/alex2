# Try to find CLHEP
#
# Once run this will define: 
# 
# CLHEP_FOUND       = system has CLHEP lib
#
# CLHEP_LIBRARIES   = full path to the libraries
#    on Unix/Linux with additional linker flags from "CLHEP-config --libs"
# 
# CMAKE_CLHEP_CXX_FLAGS  = Unix compiler flags for CLHEP, essentially "`CLHEP-config --include`"
#
# CLHEP_INCLUDE_DIR      = where to find headers 
#
# CLHEP_LINK_DIRECTORIES = link directories, useful for rpath on Unix
# CLHEP_EXE_LINKER_FLAGS = rpath on Unix
#
# JJGC modified from Felix Woelk 07/2004
# Jan Woetzel
#
# www.mip.informatik.uni-kiel.de
# --------------------------------


SET(CLHEP_CONFIG_PATH 
  $ENV{CLHEP_DIR}/bin
) 
 
FIND_PROGRAM(CLHEP_CONFIG clhep-config
  ${CLHEP_CONFIG_PATH}
)

MESSAGE("DBG CLHEP_CONFIG ${CLHEP_CONFIG}")
   
IF (CLHEP_CONFIG) 

  # ask CLHEP-config for the library dir
  # Set CLHEP_LIBRARY_DIR

  EXEC_PROGRAM(${CLHEP_CONFIG}
        ARGS --prefix
        OUTPUT_VARIABLE CLHEP_PREFIX)

  MESSAGE("DBG CLHEP_PREFIX ${CLHEP_PREFIX}")
  SET(CLHEP_INCLUDE_DIR ${CLHEP_PREFIX}/include CACHE STRING INTERNAL)
  MESSAGE("DBG CLHEP_INCLUDE_DIR ${CLHEP_INCLUDE_DIR}")
  SET(CLHEP_LIBRARY_DIR ${CLHEP_PREFIX}/lib CACHE STRING INTERNAL)
  MESSAGE("DBG CLHEP_LIBRARY_DIR ${CLHEP_LIBRARY_DIR}")

      
  # ask CLHEP-config for the library varaibles
  EXEC_PROGRAM( ${CLHEP_CONFIG}
    ARGS "--libs" 
    OUTPUT_VARIABLE CLHEP_flags )

  SET(CLHEP_LIBRARIES ${CLHEP_flags})

  MESSAGE("DBG CLHEP_LIBRARIES ${CLHEP_LIBRARIES}")


  SET(LD_LIBRARY_PATH ${LD_LIBRARY_PATH} ${CLHEP_LIBRARY_DIR})

  MESSAGE("LD_LIBRARY_PATH ${LD_LIBRARY_PATH}")

  
      
    ELSE(CLHEP_CONFIG)
      MESSAGE("FindCLHEP.cmake: CLHEP-config not found. Please set it manually. CLHEP_CONFIG=${CLHEP_CONFIG}")
    ENDIF(CLHEP_CONFIG)

IF(CLHEP_LIBRARIES)
  IF(CLHEP_INCLUDE_DIR OR CLHEP_CXX_FLAGS)

    SET(CLHEP_FOUND 1)
    
  ENDIF(CLHEP_INCLUDE_DIR OR CLHEP_CXX_FLAGS)
ENDIF(CLHEP_LIBRARIES)

