# Source: https://github.com/martinsch/pgmlink

# Modifications by Bram Custers
# b.a.custers@tue.nl

# This module finds cplex.
#
# User can give CPLEX_ROOT_DIR as a hint stored in the cmake cache.
#
# It sets the following variables:
#  CPLEX_FOUND              - Set to false, or undefined, if cplex isn't found.
#  CPLEX_INCLUDE_DIRS       - include directory
#  CPLEX_LIBRARIES          - library files


set(CPLEX_WIN_PLATFORM "")

# Try to find CPLEX directories
FIND_PATH(CPLEX_INCLUDE_DIR
        ilcplex/cplex.h
        HINTS ${CPLEX_ROOT_DIR}/cplex/include
        ${CPLEX_ROOT_DIR}/include
        PATHS ENV C_INCLUDE_PATH
        ENV C_PLUS_INCLUDE_PATH
        ENV INCLUDE_PATH
        )

MESSAGE(STATUS "INCLUDE_DIR:${CPLEX_INCLUDE_DIR}")


FIND_LIBRARY(CPLEX_LIBRARY
        NAMES cplex${CPLEX_WIN_VERSION}0 cplex
        HINTS
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_debian4.0_4.1/static_pic #unix
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_sles10_4.1/static_pic #unix
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_linux/static_pic #unix
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_osx/static_pic #osx
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_darwin/static_pic #osx
        PATHS ENV LIBRARY_PATH #unix
        ENV LD_LIBRARY_PATH #unix
        )

message("CPLEX ROOT DIR is: ${CPLEX_ROOT_DIR}")

FIND_LIBRARY(CPLEX_ILOCPLEX_LIBRARY
        ilocplex
        HINTS
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_debian4.0_4.1/static_pic #unix
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_sles10_4.1/static_pic #unix
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_linux/static_pic #unix
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_osx/static_pic #osx
        ${CPLEX_ROOT_DIR}/cplex/lib/x86-64_darwin/static_pic #osx
        PATHS ENV LIBRARY_PATH
        ENV LD_LIBRARY_PATH
        )

message(STATUS "ILOCPLEX Library: ${CPLEX_ILOCPLEX_LIBRARY}")

FIND_PATH(CPLEX_BIN_DIR
        cplex
        HINTS ${CPLEX_ROOT_DIR}/cplex/bin/x86-64_sles10_4.1 #unix
        ${CPLEX_ROOT_DIR}/cplex/bin/x86-64_debian4.0_4.1 #unix
        ${CPLEX_ROOT_DIR}/cplex/bin/x86-64_linux #unix
        ${CPLEX_ROOT_DIR}/cplex/bin/x86-64_osx #osx
        ${CPLEX_ROOT_DIR}/cplex/bin/x86-64_darwin #osx
        ENV LIBRARY_PATH
        ENV LD_LIBRARY_PATH
        )

message(STATUS "CPLEX Bin Dir: ${CPLEX_BIN_DIR}")

# Apply the standard find package function.
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(CPLEX DEFAULT_MSG
        CPLEX_LIBRARY CPLEX_INCLUDE_DIR CPLEX_ILOCPLEX_LIBRARY)


# Set the output variables for CPLEX
IF(CPLEX_FOUND)
    SET(CPLEX_INCLUDE_DIRS ${CPLEX_INCLUDE_DIR})
    SET(CPLEX_LIBRARIES ${CPLEX_ILOCPLEX_LIBRARY} ${CPLEX_LIBRARY} )
    IF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
        SET(CPLEX_LIBRARIES "${CPLEX_LIBRARIES};m;pthread")
    ENDIF(CMAKE_SYSTEM_NAME STREQUAL "Linux")
else()
    message(FATAL_ERROR "CPLEX not found. Please install it; if installed, try defining CPLEX_ROOT_DIR.")
ENDIF(CPLEX_FOUND)

# Mark variable advanced, but inspectable
MARK_AS_ADVANCED(CPLEX_LIBRARY CPLEX_INCLUDE_DIR CPLEX_ILOCPLEX_LIBRARY)

# Create an interface library that will contain the relevant data to link against.
# This way, you don't have to do all the target_* yourself for directories and libraries.
add_library(Cplex INTERFACE IMPORTED)

# Add the Cplex directories
set_property(TARGET Cplex APPEND PROPERTY INTERFACE_INCLUDE_DIRECTORIES
        ${CPLEX_INCLUDE_DIRS}
        )

# Add the Cplex libraries
set_property(TARGET Cplex APPEND PROPERTY INTERFACE_LINK_LIBRARIES
        ${CPLEX_LIBRARIES} # Generator expression for this: SET(CPLEX_LIBRARIES "${CPLEX_LIBRARIES};m;pthread")
        )

# Add Dll copy target for windows to copy CPLEX dll's to the output folder (in case CPLEX is not on the system path.)
if(WIN32)
    add_custom_target(CplexDllCopy
            # todo: check if debug and release folder exist
            # debug version
            COMMAND ${CMAKE_COMMAND} -E copy ${CPLEX_BIN_DIR}/cplex${CPLEX_WIN_VERSION}0.dll          ${CMAKE_BINARY_DIR}
            )
    add_dependencies(Cplex CplexDllCopy)
endif()
