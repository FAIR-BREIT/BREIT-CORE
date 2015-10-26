 ################################################################################
 #    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    #
 #                                                                              #
 #              This software is distributed under the terms of the             # 
 #         GNU Lesser General Public Licence version 3 (LGPL) version 3,        #  
 #                  copied verbatim in the file "LICENSE"                       #
 ################################################################################
# Create a library called "libbreit_utils" 

set(INCLUDE_DIRECTORIES ${CMAKE_SOURCE_DIR}/breit-utils/src ${CMAKE_SOURCE_DIR}/breit-utils/run)

set(SYSTEM_INCLUDE_DIRECTORIES ${SYSTEM_INCLUDE_DIRECTORIES})

include_directories(${INCLUDE_DIRECTORIES})
include_directories(${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES ${Boost_LIBRARY_DIRS}) 

link_directories(${LINK_DIRECTORIES})

set(SRCS
      src/options_manager.cxx
      src/logger.cxx
    )

set(LIBRARY_NAME breit_utils)
set(DEPENDENCIES 
    boost_log 
    boost_log_setup 
    boost_thread 
    boost_date_time 
    boost_filesystem 
    boost_system
    boost_program_options
    pthread
)

GENERATE_LIBRARY()

set(EXE_NAME runtestLogger)
set(SRCS run/testLogger.cxx)
set(DEPENDENCIES breit_utils)
GENERATE_EXECUTABLE()
