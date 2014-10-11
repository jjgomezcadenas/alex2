# -*- cmake -*-

FILE(GLOB algos_src *.cpp *.cxx )

# # Finally, create it!
ADD_LIBRARY(alexAlgos SHARED ${algos_src})
SET_TARGET_PROPERTIES(alexAlgos PROPERTIES VERSION 1.1.0 SOVERSION 1)
TARGET_LINK_LIBRARIES(alexAlgos alex ${LOG4CPP_LIBRARIES} ${IRENE_LIBRARIES} ${PAOLINA_LIBRARIES} 
	${RECPACK_LIBRARIES} ${ROOT_LIBRARIES} ${CLHEP_LIBRARIES} ${TINYXML2_LIBRARIES})

SET(CMAKE_VERBOSE_MAKEFILE ON)

# Installation destinations (just for the lib, not the python module)
INSTALL(TARGETS alexAlgos LIBRARY DESTINATION lib)
