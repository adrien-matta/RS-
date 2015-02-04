cmake_minimum_required (VERSION 2.6)
project (RS++)
 include("FindROOT.cmake") 
# The version number.
set (RS_VERSION_MAJOR 1)
set (RS_VERSION_MINOR 0)

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_FLAGS "-std=c++11 -stdlib=libc++")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
set(CMAKE_INCLUDE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/include )

#link_directories( ${LINK_DIRECTORIES})
# add the binary tree to the search path for include files
# so that we will find RSConfig.h
#include_directories(${CMAKE_INCLUDE_OUTPUT_DIRECTORY})

set(INCLUDE_DIRECTORIES ${ROOT_INCLUDE_DIR} )
include_directories( ${INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES ${ROOT_LIBRARY_DIR})
link_directories( ${LINK_DIRECTORIES})

# add lib
set(RS_SUB_DIRECTORIES ShellModel Reaction)

# Macro that add the sub dir and the include file install path
foreach(subdir ${RS_SUB_DIRECTORIES})
 add_subdirectory (${subdir})  
 
file(GLOB incfiles "${subdir}/*.h")
  foreach(incfile ${incfiles})
  install (FILES ${incfile}
           DESTINATION ${CMAKE_BINARY_DIR}/include)
  endforeach(incfile)
endforeach(subdir)
