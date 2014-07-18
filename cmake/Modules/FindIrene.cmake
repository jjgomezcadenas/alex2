# Try to find irene
#
# Once run this will define: 
# 
# IRENE_FOUND       = system has IRENE lib
#
# IRENE_LIBRARIES   = full path to the libraries
#    on Unix/Linux with additional linker flags from "irene-config --libs"
# 
# CMAKE_IRENE_CXX_FLAGS  = Unix compiler flags for IRENE, essentially "`irene-config --include`"
#
# IRENE_INCLUDE_DIR      = where to find headers 
#
# IRENE_LINK_DIRECTORIES = link directories, useful for rpath on Unix
# IRENE_EXE_LINKER_FLAGS = rpath on Unix
#
# JJGC modified from Felix Woelk 07/2004
# Jan Woetzel
#
# www.mip.informatik.uni-kiel.de
# --------------------------------


SET(IRENE_CONFIG_PATH 
  $ENV{IRENE}/bin
) 
 
FIND_PROGRAM(IRENE_CONFIG irene-config
  ${IRENE_CONFIG_PATH}
)

MESSAGE("DBG IRENE_CONFIG ${IRENE_CONFIG}")
   
IF (IRENE_CONFIG) 

  # ask irene-config for the library dir
  # Set IRENE_LIBRARY_DIR

  EXEC_PROGRAM( ${IRENE_CONFIG}
    ARGS "--libdir"
    OUTPUT_VARIABLE IRENE_LIBRARY_DIR)
    
  MESSAGE("DBG IRENE_LIBRARY_DIR ${IRENE_LIBRARY_DIR}")

  # ask irene-config for the include dir
  #EXEC_PROGRAM( ${IRENE_CONFIG}
  #  ARGS "--include" 
  #  OUTPUT_VARIABLE irene_headers )
  #SET(IRENE_INCLUDE_DIR ${irene_headers})

  # set INCLUDE_DIRS to prefix+include
      EXEC_PROGRAM(${IRENE_CONFIG}
        ARGS --prefix
        OUTPUT_VARIABLE IRENE_PREFIX)
      SET(IRENE_INCLUDE_DIR ${IRENE_PREFIX}/include CACHE STRING INTERNAL)


  MESSAGE("DBG IRENE_INCLUDE_DIR ${IRENE_INCLUDE_DIR}")
      
  # ask irene-config for the library varaibles
  EXEC_PROGRAM( ${IRENE_CONFIG}
    ARGS "--libs" 
    OUTPUT_VARIABLE irene_flags )

  SET(IRENE_LIBRARIES ${irene_flags})

  MESSAGE("DBG IRENE_LIBRARIES ${IRENE_LIBRARIES}")

  # Set IRENE_INCLUDES
  SET( IRENE_INCLUDES ${IRENE_INCLUDE_DIR})

  MESSAGE("DBG IRENE_INCLUDES ${IRENE_INCLUDES}")

  SET(LD_LIBRARY_PATH ${LD_LIBRARY_PATH} ${IRENE_LIBRARY_DIR})

  MESSAGE("LD_LIBRARY_PATH ${LD_LIBRARY_PATH}")

  #SET(IRENE_CXX_FLAGS "`${IRENE_CONFIG} --include`")

  #MESSAGE("DBG IRENE_CXX_FLAGS ${IRENE_CXX_FLAGS}")
      
  # set INCLUDE_DIRS to prefix+include
  #EXEC_PROGRAM(${IRENE_CONFIG}
  #      ARGS --prefix
  #      OUTPUT_VARIABLE IRENE_PREFIX)
  #      SET(IRENE_INCLUDE_DIR ${IRENE_PREFIX}/include CACHE STRING INTERNAL)


  #MESSAGE("IRENE_LIBRARY_DIR ${IRENE_LIBRARY_DIR}")

  # set link libraries and link flags
   # SET(IRENE_LIBRARIES "`${IRENE_CONFIG} --libs`")
      
      # extract link dirs for rpath  
  #    EXEC_PROGRAM(${IRENE_CONFIG}
  #      ARGS --libdir
   #     OUTPUT_VARIABLE IRENE_CONFIG_LIBS )

      # split off the link dirs (for rpath)
      # use regular expression to match wildcard equivalent "-L*<endchar>"
      # with <endchar> is a space or a semicolon
  #    STRING(REGEX MATCHALL "[-][L]([^ ;])+" 
  #      IRENE_LINK_DIRECTORIES_WITH_PREFIX 
  #      "${IRENE_CONFIG_LIBS}" )
      #      MESSAGE("DBG  IRENE_LINK_DIRECTORIES_WITH_PREFIX=${IRENE_LINK_DIRECTORIES_WITH_PREFIX}")

      # remove prefix -L because we need the pure directory for LINK_DIRECTORIES
      
  #    IF (IRENE_LINK_DIRECTORIES_WITH_PREFIX)
  #      STRING(REGEX REPLACE "[-][L]" "" IRENE_LINK_DIRECTORIES ${IRENE_LINK_DIRECTORIES_WITH_PREFIX} )
  #    ENDIF (IRENE_LINK_DIRECTORIES_WITH_PREFIX)
  #    SET(IRENE_EXE_LINKER_FLAGS "-Wl,-rpath,${IRENE_LINK_DIRECTORIES}" CACHE STRING INTERNAL)
      #      MESSAGE("DBG  IRENE_LINK_DIRECTORIES=${IRENE_LINK_DIRECTORIES}")
      #      MESSAGE("DBG  IRENE_EXE_LINKER_FLAGS=${IRENE_EXE_LINKER_FLAGS}")

      #      ADD_DEFINITIONS("-DHAVE_IRENE")
      #      SET(IRENE_DEFINITIONS "-DHAVE_IRENE")
  #    MARK_AS_ADVANCED(
  #      IRENE_CXX_FLAGS
  #      IRENE_INCLUDE_DIR
  #      IRENE_LIBRARIES
  #      IRENE_LINK_DIRECTORIES
  #      IRENE_DEFINITIONS
	#)
  #    MESSAGE(STATUS "Using IRENE from ${IRENE_PREFIX}")
      
    ELSE(IRENE_CONFIG)
      MESSAGE("Findirene.cmake: irene-config not found. Please set it manually. IRENE_CONFIG=${IRENE_CONFIG}")
    ENDIF(IRENE_CONFIG)

IF(IRENE_LIBRARIES)
  IF(IRENE_INCLUDE_DIR OR IRENE_CXX_FLAGS)

    SET(IRENE_FOUND 1)
    
  ENDIF(IRENE_INCLUDE_DIR OR IRENE_CXX_FLAGS)
ENDIF(IRENE_LIBRARIES)

