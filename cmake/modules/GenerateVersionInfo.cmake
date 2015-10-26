 ################################################################################
 #    Copyright (C) 2014 GSI Helmholtzzentrum fuer Schwerionenforschung GmbH    #
 #                                                                              #
 #              This software is distributed under the terms of the             # 
 #         GNU Lesser General Public Licence version 3 (LGPL) version 3,        #  
 #                  copied verbatim in the file "LICENSE"                       #
 ################################################################################
 # Configure FairVersion.h
 # ------------------------------

 Find_Package(Git)

 If(GIT_FOUND AND EXISTS "${SOURCE_DIR}/.git")
   Execute_Process(COMMAND ${GIT_EXECUTABLE} describe --always
                   OUTPUT_VARIABLE BREIT_GIT_VERSION
                   OUTPUT_STRIP_TRAILING_WHITESPACE
                   WORKING_DIRECTORY ${SOURCE_DIR}
                  )
   Execute_Process(COMMAND ${GIT_EXECUTABLE} log -1 --format=%cd
                   OUTPUT_VARIABLE BREIT_GIT_DATE
                   OUTPUT_STRIP_TRAILING_WHITESPACE
                   WORKING_DIRECTORY ${SOURCE_DIR}
                  )
   Message(STATUS "BREIT Version - ${BREIT_GIT_VERSION} from - ${BREIT_GIT_DATE}")
   if(BREIT)
     Configure_File(${BREIT}/scripts/FairVersion.h.tmp ${BINARY_DIR}/FairVersion.h @ONLY)
   else(BREIT)  
     Configure_File(${SOURCE_DIR}/cmake/scripts/FairVersion.h.tmp ${BINARY_DIR}/FairVersion.h @ONLY)
   endif(BREIT)
  
 Else()
   if(BREIT)
     Configure_File(${BREIT}/scripts/FairVersion.h.default ${BINARY_DIR}/FairVersion.h COPYONLY)
   else(BREIT) 
     Configure_File(${SOURCE_DIR}/cmake/scripts/FairVersion.h.default ${BINARY_DIR}/FairVersion.h COPYONLY)
   endif(BREIT)
 EndIf()

