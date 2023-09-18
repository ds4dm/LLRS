cmake_minimum_required(VERSION 3.18)
project(Assortment_Optimization)

set(Rossobianco true)
message(STATUS "Rossobianco variable=${Rossobianco}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ldl")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -ldl")
######## Cplex Dependency
if (Rossobianco)
    set(CPLEX_ROOT_DIR "/home/ibm/cplex-studio/20.1")
elseif(Cedar)
    set(CPLEX_ROOT_DIR "/home/clsol/opt/ibm/ILOG/CPLEX_Studio201")
else()
    set(CPLEX_ROOT_DIR "/Users/claudiosole/Applications/IBM/ILOG/CPLEX_Studio1271")
endif()
find_package(CPLEX REQUIRED)
#include_directories(${CPLEX_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_C_STANDARD 11)

include_directories(AssortOpt/AssortOpt ${CPLEX_INCLUDE_DIRS})
add_executable(ASOpt AssortOpt/AssortOpt/Source.cpp)
add_executable(ASOpt_MMNL AssortOpt/AssortOpt/Source_MMNL.cpp)
target_link_libraries(ASOpt ${CPLEX_LIBRARIES})
target_link_libraries(ASOpt_MMNL ${CPLEX_LIBRARIES})
