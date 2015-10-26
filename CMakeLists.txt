# The name of our project is “BREIT”.  CMakeLists files in this project can
# refer to the root source directory of the project as ${BREIT_SOURCE_DIR} 
# or as ${CMAKE_SOURCE_DIR} and to the root binary directory of the project as 
# ${BREIT_BINARY_DIR} or ${CMAKE_BINARY_DIR}.
# This difference is important for the base classes which are in BREIT 
# and the experiment part.

# Check if cmake has the required version
CMAKE_MINIMUM_REQUIRED(VERSION  2.8.2 FATAL_ERROR)



### CMP0025   Compiler id for Apple Clang is now AppleClang.
### CMP0042   MACOSX_RPATH is enabled by default.

foreach(p
  CMP0025 # CMake 3.0
  CMP0042 # CMake 3.0
  CMP0054
  )
  if(POLICY ${p})
  cmake_policy(SET ${p} NEW)
  endif()
endforeach()


# Set name of our project to “BREIT". Has to be done
# after check of cmake version since this is a new feature                  
project(BREIT)

#In case you need Fortran 
#ENABLE_LANGUAGE(Fortran)



# where to look first for cmake modules, before ${CMAKE_ROOT}/Modules/ 
# is checked
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

Set(CheckSrcDir "${CMAKE_SOURCE_DIR}/cmake/checks")
include(CheckCXX11Features)


# Load some basic macros which are needed later on
include(FairMacros)
include(WriteConfigFile)
include(Dart)
include(CheckCompiler)

#Check the compiler and set the compile and link flags

#Check the compiler and set the compile and link flags
Check_Compiler()

set(LIBRARY_OUTPUT_PATH "${CMAKE_BINARY_DIR}/lib")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_BINARY_DIR}/bin")

Option(USE_PATH_INFO "Information from PATH and LD_LIBRARY_PATH are
used." OFF)
If(USE_PATH_INFO)
  Set(PATH $ENV{PATH})
  If (APPLE)
    Set(LD_LIBRARY_PATH $ENV{DYLD_LIBRARY_PATH})
  Else (APPLE)
    Set(LD_LIBRARY_PATH $ENV{LD_LIBRARY_PATH})
  EndIf (APPLE)
EndIf(USE_PATH_INFO)

# Check if the user wants to build the project in the source
# directory and if the install directory is different from the build
# directory
CHECK_OUT_OF_SOURCE_BUILD()
CHECK_INSTALL_DIRECTORY()

# Check if we are on an UNIX system. If not stop with an error
# message
IF(NOT UNIX)
  MESSAGE(FATAL_ERROR "You're not on an UNIX system. The project was up to now only tested on UNIX systems, so we break here. If you want to go on please edit the CMakeLists.txt in the source directory.")
ENDIF(NOT UNIX)  


#Check if necessary environment variables are set
#If not stop execution unless modular build is activated
Option(BREIT_MODULAR_BUILD "Modular build without environment variables" OFF)
IF(NOT BREIT_MODULAR_BUILD)
  IF(NOT DEFINED ENV{SIMPATH})
    MESSAGE(FATAL_ERROR "You did not define the environment variable SIMPATH which is nedded to find the external packages. Please set this variable and execute cmake again.")
  ENDIF(NOT DEFINED ENV{SIMPATH})
  STRING(REGEX MATCHALL "[^:]+" PATH $ENV{PATH})
  SET(SIMPATH $ENV{SIMPATH})
endif(NOT BREIT_MODULAR_BUILD)

# Check if the external packages are installed into a separate install
# directory
CHECK_EXTERNAL_PACKAGE_INSTALL_DIR()

# searches for needed packages
# REQUIRED means that cmake will stop if this packages are not found
# For example the framework can run without GEANT4, but ROOT is
# mandatory

find_package(ROOT 5.32.00)


If(FAIRSOFT_EXTERN)
  set(BOOST_ROOT ${SIMPATH})
  set(GSL_DIR ${SIMPATH})
Else(FAIRSOFT_EXTERN)
  set(BOOST_ROOT ${SIMPATH}/basics/boost)
  set(GSL_DIR ${SIMPATH}/basics/gsl)
EndIf(FAIRSOFT_EXTERN)

Set(Boost_NO_SYSTEM_PATHS TRUE)
Set(Boost_NO_BOOST_CMAKE TRUE)
Message("-- Looking for Boost ...")
# If an older version of boost is found both of the variables below are
# cached and in a second cmake run, a good boost version is found even 
# if the version is to old. 
# To overcome this problem both variables are cleared before checking
# for boost.
Unset(Boost_INCLUDE_DIR CACHE)
Unset(Boost_LIBRARY_DIRS CACHE)
find_package(Boost 1.41 COMPONENTS thread system timer program_options log)
If (Boost_FOUND)
  Set(Boost_Avail 1)
  Set(LD_LIBRARY_PATH ${LD_LIBRARY_PATH} ${Boost_LIBRARY_DIR})
Else (Boost_FOUND)
  Set(Boost_Avail 0)
EndIf (Boost_FOUND)


find_package(LAPACK)
#find_package(BoostNumBinding)

# Set the library version in the main CMakeLists.txt
SET(BREIT_MAJOR_VERSION 0)
SET(BREIT_MINOR_VERSION 0)
SET(BREIT_PATCH_VERSION 0)
SET(BREIT_VERSION "${BREIT_MAJOR_VERSION}.${BREIT_MINOR_VERSION}.${BREIT_PATCH_VERSION}")
SET(BREIT_LIBRARY_PROPERTIES ${BREIT_LIBRARY_PROPERTIES}
    VERSION "${BREIT_VERSION}"
    SOVERSION "${BREIT_MAJOR_VERSION}"
    SUFFIX ".so"
)

# Load the macro to generate the needed version info at compile time
Generate_Version_Info()

# Check if the compiler support specific C++11 features
# Up to now this is only a check since the code does not use
# any of the features of the new standard
include(CheckCXX11Features)

IF(HAS_CXX11_SHAREDPOINTER)
  Add_Definitions(-DHAS_SHAREDPOINTER)
ENDIF(HAS_CXX11_SHAREDPOINTER)

# Recurse into the given subdirectories.  This does not actually
# cause another cmake executable to run.  The same process will walk through
# the project's entire directory structure.
#SetBasicVariables()

Set(SYSTEM_INCLUDE_DIRECTORIES
    ${Boost_INCLUDE_DIRS}
  )  

#add_subdirectory(breit-doc)

If (Boost_FOUND)
    add_subdirectory(breit-utils)
    add_subdirectory(breit-tests)
endif(Boost_FOUND)

WRITE_CONFIG_FILE(config.sh)
 
configure_file(${CMAKE_SOURCE_DIR}/CTestCustom.cmake 
               ${CMAKE_BINARY_DIR}/CTestCustom.cmake
              )


