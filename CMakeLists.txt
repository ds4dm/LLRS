cmake_minimum_required(VERSION 3.18)
project(Assortment_Optimization)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ldl")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -ldl")
######## Cplex Dependency
set(CPLEX_ROOT_DIR "...") # Add Cplex dependency here
find_package(CPLEX REQUIRED)
include_directories(${CPLEX_INCLUDE_DIRS})
IF(APPLE)
    SET(CMAKE_OSX_ARCHITECTURES "x86_64" CACHE STRING "Build architectures for Mac OS X" FORCE)
ENDIF(APPLE)
set(CMAKE_CXX_STANDARD 17)
include_directories(AssortOpt/AssortOpt ${CPLEX_INCLUDE_DIRS})
add_executable(ASOpt AssortOpt/AssortOpt/Source.cpp)
add_executable(ASOpt_MMNL AssortOpt/AssortOpt/Source_MMNL.cpp)
target_link_libraries(ASOpt ${CPLEX_LIBRARIES})
target_link_libraries(ASOpt_MMNL ${CPLEX_LIBRARIES})




