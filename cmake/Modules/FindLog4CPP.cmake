# Try to find LOG4CPP
#
# Once run this will define: 
# 
# LOG4CPP_FOUND       = system has LOG4CPP lib
#
# LOG4CPP_LIBRARIES   = full path to the libraries
#    on Unix/Linux with additional linker flags from "LOG4CPP-config --libs"
# 
# CMAKE_LOG4CPP_CXX_FLAGS  = Unix compiler flags for LOG4CPP, essentially "`LOG4CPP-config --include`"
#
# LOG4CPP_INCLUDE_DIR      = where to find headers 
#
# LOG4CPP_LINK_DIRECTORIES = link directories, useful for rpath on Unix
# LOG4CPP_EXE_LINKER_FLAGS = rpath on Unix
#
# JJGC modified from Felix Woelk 07/2004
# Jan Woetzel
#
# www.mip.informatik.uni-kiel.de
# --------------------------------


SET(LOG4CPP_CONFIG_PATH 
  $ENV{LOG4CPP_DIR}/bin
) 
 
FIND_PROGRAM(LOG4CPP_CONFIG log4cpp-config
  ${LOG4CPP_CONFIG_PATH}
)

MESSAGE("DBG LOG4CPP_CONFIG ${LOG4CPP_CONFIG}")
   
IF (LOG4CPP_CONFIG) 

  # ask LOG4CPP-config for the library dir
  # Set LOG4CPP_LIBRARY_DIR

  EXEC_PROGRAM( ${LOG4CPP_CONFIG}
    ARGS "--libdir"
    OUTPUT_VARIABLE LOG4CPP_LIBRARY_DIR)
    
  MESSAGE("DBG LOG4CPP_LIBRARY_DIR ${LOG4CPP_LIBRARY_DIR}")

  # ask LOG4CPP-config for the include dir
  #EXEC_PROGRAM( ${LOG4CPP_CONFIG}
  #  ARGS "--cflags" 
  #  OUTPUT_VARIABLE LOG4CPP_headers )
  #SET(LOG4CPP_INCLUDE_DIR ${LOG4CPP_headers})

  # set INCLUDE_DIRS to prefix+include
      EXEC_PROGRAM(${LOG4CPP_CONFIG}
        ARGS --prefix
        OUTPUT_VARIABLE LOG4CPP_PREFIX)
      SET(LOG4CPP_INCLUDE_DIR ${LOG4CPP_PREFIX}/include CACHE STRING INTERNAL)

  MESSAGE("DBG LOG4CPP_INCLUDE_DIR ${LOG4CPP_INCLUDE_DIR}")
      
  # ask log4cpp-config for the library varaibles
  EXEC_PROGRAM( ${LOG4CPP_CONFIG}
    ARGS "--libs" 
    OUTPUT_VARIABLE LOG4CPP_flags )

  SET(LOG4CPP_LIBRARIES ${LOG4CPP_flags})

  MESSAGE("DBG LOG4CPP_LIBRARIES ${LOG4CPP_LIBRARIES}")

  # Set LOG4CPP_INCLUDES
  SET( LOG4CPP_INCLUDES ${LOG4CPP_INCLUDE_DIR})

  MESSAGE("DBG LOG4CPP_INCLUDES ${LOG4CPP_INCLUDES}")

  SET(LD_LIBRARY_PATH ${LD_LIBRARY_PATH} ${LOG4CPP_LIBRARY_DIR})

  MESSAGE("LD_LIBRARY_PATH ${LD_LIBRARY_PATH}")

  
    ELSE(LOG4CPP_CONFIG)
      MESSAGE("FindLOG4CPP.cmake: LOG4CPP-config not found. Please set it manually. LOG4CPP_CONFIG=${LOG4CPP_CONFIG}")
    ENDIF(LOG4CPP_CONFIG)

IF(LOG4CPP_LIBRARIES)
  IF(LOG4CPP_INCLUDE_DIR OR LOG4CPP_CXX_FLAGS)

    SET(LOG4CPP_FOUND 1)
    
  ENDIF(LOG4CPP_INCLUDE_DIR OR LOG4CPP_CXX_FLAGS)
ENDIF(LOG4CPP_LIBRARIES)

