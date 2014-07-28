# -*- cmake -*-

FILE(GLOB example_src *.cpp *.cxx )

# # Finally, create it!
ADD_LIBRARY(alexExample SHARED ${example_src})
SET_TARGET_PROPERTIES(alexExample PROPERTIES VERSION 1.1.0 SOVERSION 1)
TARGET_LINK_LIBRARIES(alexExample alex ${LOG4CPP_LIBRARIES} ${IRENE_LIBRARIES} ${PAOLINA_LIBRARIES} 
	${RECPACK_LIBRARIES} ${ROOT_LIBRARIES} ${CLHEP_LIBRARIES} ${TINYXML2_LIBRARIES})

SET(CMAKE_VERBOSE_MAKEFILE ON)

# Installation destinations (just for the lib, not the python module)
INSTALL(TARGETS alexExample LIBRARY DESTINATION lib)
