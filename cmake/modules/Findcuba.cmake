 ################################################################################
 #    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    #
 #                                                                              #
 #              This software is distributed under the terms of the             # 
 #         GNU Lesser General Public Licence version 3 (LGPL) version 3,        #  
 #                  copied verbatim in the file "LICENSE"                       #
 ################################################################################
# - Try to find BAT
# Once found, this will define
#
#  CUBA_FOUND - system has CUBA
#  CUBA_INCLUDE_DIR - the CUBA include directories
#  CUBA_LIBRARY_DIR - the CUBA lib directories
#  CUBA_LIBRARIES_STATIC 
message(STATUS "Looking for CUBA...")

set(_HEADERS cuba.h)

set(SEARCH_DIRS /usr/local/lib /usr/local/include /opt/local/include /opt/local/lib)

if (APPLE)
  set(LIBCUBA_SHARED libcuba.dylib)
else (APPLE)
  set(LIBCUBA_SHARED libcuba.so)
endIf (APPLE)

set(LIBCUBA_STATIC libcuba.a)

find_path(CUBA_INCLUDE_DIR NAMES ${_HEADERS} 
  PATHS ${SEARCH_DIRS}
  DOC   "Path to CUBA include header files."
)

find_path(CUBA_LIBRARY_DIR NAMES libcuba.a
  PATHS ${SEARCH_DIRS}
  DOC   "Path to CUBA lib files."
)

find_library(CUBA_LIBRARY_SHARED NAMES libcuba.dyld
  PATHS ${CUBA_LIBRARY_DIR}
  DOC   "Path to ${LIBCUBA_SHARED}."
)

find_library(CUBA_LIBRARY_STATIC NAMES libcuba.a
  PATHS ${CUBA_LIBRARY_DIR}
  NO_DEFAULT_PATH
  DOC   "Path to ${LIBCUBA_STATIC}."
)

#if(CUBA_INCLUDE_DIR AND CUBA_LIBRARY_SHARED AND CUBA_LIBRARY_STATIC)
#  set(CUBA_FOUND true)
#else(CUBA_INCLUDE_DIR AND CUBA_LIBRARY_SHARED AND CUBA_LIBRARY_STATIC)
#  set(CUBA_FOUND false)
#endif(CUBA_INCLUDE_DIR AND CUBA_LIBRARY_SHARED AND CUBA_LIBRARY_STATIC)


if(CUBA_INCLUDE_DIR AND CUBA_LIBRARY_STATIC)
  set(CUBA_FOUND true)
else(CUBA_INCLUDE_DIR AND CUBA_LIBRARY_STATIC)
  set(CUBA_FOUND false)
endif(CUBA_INCLUDE_DIR AND CUBA_LIBRARY_STATIC)


if(CUBA_FOUND)
  add_definitions(-DCUBA)
  set(CUBA_LIBRARIES "${CUBA_LIBRARY_STATIC};${CUBA_LIBRARY_SHARED}")
  if(NOT CUBA_FIND_QUIETLY)
    #message(STATUS "Looking for CUBA... - Found ${CUBA_LIBRARIES}")
    message(STATUS "Looking for CUBA... - Found include dir ${CUBA_INCLUDE_DIR}")
    message(STATUS "Looking for CUBA... - Found library dir ${CUBA_LIBRARY_DIR}")
  endif(NOT CUBA_FIND_QUIETLY)
else(CUBA_FOUND)
  if(NOT CUBA_FIND_QUIETLY)
    if(CUBA_FIND_REQUIRED)
      message(FATAL_ERROR "Looking for CUBA... - Not found")
    else(CUBA_FIND_REQUIRED)
      message(STATUS "Looking for CUBA... - Not found")
    endif(CUBA_FIND_REQUIRED)
  endif(NOT CUBA_FIND_QUIETLY)
endif(CUBA_FOUND)

#mark_as_advanced(CUBA_INCLUDE_DIR CUBA_LIBRARY_DIR CUBA_LIBRARY_SHARED CUBA_LIBRARY_STATIC)
mark_as_advanced(CUBA_INCLUDE_DIR CUBA_LIBRARY_DIR CUBA_LIBRARY_STATIC)
