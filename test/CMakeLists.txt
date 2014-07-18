CMAKE_MINIMUM_REQUIRED(VERSION 2.8)
# Make sure the compiler can find include files from our  library. 
include_directories (${ALEX_SOURCE_DIR}) 

# Make sure the linker can find the  library once it is built. 
link_directories (${ALEX_SOURCE_DIR}) 

#add_executable (testAConfig testAconfig.cpp) 
#add_executable (testISvc testISvc.cpp) 
#add_executable (testPaolina testPaolina.cpp) 
add_executable (testXml testXml.cpp)
add_executable (testAlex testAlex.cpp) 



# Link the executable to the  library. 
#target_link_libraries (testAConfig alex)
#target_link_libraries (testISvc alex)
#target_link_libraries (testPaolina alex)
target_link_libraries (testXml alex)
target_link_libraries (testAlex alex)
