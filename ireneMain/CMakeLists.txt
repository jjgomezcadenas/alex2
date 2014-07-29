CMAKE_MINIMUM_REQUIRED(VERSION 2.8)


# Make sure the compiler can find include files from our  library. 
include_directories (${ALEX_SOURCE_DIR}) 

# Make sure the linker can find the  library once it is built. 
link_directories (${ALEX_SOURCE_DIR} ${ALEX_SOURCE_DIR}/example) 

add_executable (alexMainExample alexMainExample.cpp)
target_link_libraries (alexMainExample alex alexExample)


