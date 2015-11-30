

set(SCRIPT_CMD1 "\${filename##*.}")
set(SCRIPT_CMD2 "\${filename%.*}")

configure_file( ${CMAKE_SOURCE_DIR}/breit-exec/macro/breit-eq.sh.in ${CMAKE_BINARY_DIR}/bin/breit-eq.sh )
configure_file( ${CMAKE_SOURCE_DIR}/breit-exec/macro/breit-plot.sh.in ${CMAKE_BINARY_DIR}/bin/breit-plot.sh )
configure_file( ${CMAKE_SOURCE_DIR}/breit-exec/macro/breit-web.sh.in ${CMAKE_BINARY_DIR}/bin/breit-web.sh )

configure_file( ${CMAKE_SOURCE_DIR}/breit-exec/config/breit-config-save-to-file.cfg.in ${CMAKE_SOURCE_DIR}/data/config/breit-config-save-to-file.cfg )

configure_file( ${CMAKE_SOURCE_DIR}/breit-exec/config/breit-config-no-save-to-file.cfg.in ${CMAKE_SOURCE_DIR}/data/config/breit-config-no-save-to-file.cfg )

configure_file( ${CMAKE_SOURCE_DIR}/breit-exec/config/breit-config-webbreit.cfg.in ${CMAKE_SOURCE_DIR}/data/config/breit-config-webbreit.cfg )

EXEC_PROGRAM(/bin/chmod ARGS "u+x  ${CMAKE_BINARY_DIR}/bin/breit-eq.sh")
EXEC_PROGRAM(/bin/chmod ARGS "u+x  ${CMAKE_BINARY_DIR}/bin/breit-plot.sh")
EXEC_PROGRAM(/bin/chmod ARGS "u+x  ${CMAKE_BINARY_DIR}/bin/breit-web.sh")


set(INCLUDE_DIRECTORIES
${CMAKE_SOURCE_DIR}/breit-utils/src
${CMAKE_SOURCE_DIR}/breit-equations/src
${CMAKE_SOURCE_DIR}/breit-exec/policy-impl
)

Set(SYSTEM_INCLUDE_DIRECTORIES 
  ${Boost_INCLUDE_DIRS}
  ${CMAKE_SOURCE_DIR}/breit-matrix-operations
  ${CMAKE_SOURCE_DIR}/breit-matrix-operations/src
  )

if(LAPACK_FOUND)
  set(INCLUDE_DIRECTORIES 
      #${BNB_INCLUDE_DIR}
      ${INCLUDE_DIRECTORIES}
     )
  if(ROOT_FOUND)
        set(INCLUDE_DIRECTORIES 
            ${ROOT_INCLUDE_DIR}
            ${INCLUDE_DIRECTORIES}
           )
  endif(ROOT_FOUND)
endif(LAPACK_FOUND)

include_directories(${INCLUDE_DIRECTORIES})
include_directories(SYSTEM ${SYSTEM_INCLUDE_DIRECTORIES})

set(LINK_DIRECTORIES ${Boost_LIBRARY_DIRS})
set(LINK_DIRECTORIES ${ROOT_LIBRARY_DIR} ${LINK_DIRECTORIES}) 

link_directories(${LINK_DIRECTORIES})


Set(EXE_NAME runSolveSteadyEq)
Set(SRCS run/runSolveSystemAtEquilibrium.cxx)
Set(DEPENDENCIES breit_utils blas lapack gfortran)
GENERATE_EXECUTABLE()

if(LAPACK_FOUND)
  Set(EXE_NAME runSolveSteadyEqLapack)
  Set(SRCS 
    run/runSolveSystemAtEquilibriumLapack.cxx
  )
  #set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
  #set(LIBS_MAP_LAPACK "LAPACK")
  Set(DEPENDENCIES breit_utils blas lapack gfortran ${LAPACK_LIBRARIES})
  GENERATE_EXECUTABLE()


  ## ROOT GUI
  if(ROOT_FOUND)
  
      #set(LINK_DIRECTORIES ${ROOT_LIBRARY_DIR} ${LINK_DIRECTORIES}) 

      #link_directories(${LINK_DIRECTORIES})

      Set(SRCS ${CMAKE_SOURCE_DIR}/breit-exec/policy-impl/handle_root_signal.cxx)

      Set(LINKDEF ${CMAKE_SOURCE_DIR}/breit-exec/policy-impl/sig_handlerLinkDef.h)
      Set(LIBRARY_NAME handle_root_signal)

      if (ROOT_FOUND_VERSION LESS 59999)
        Set(DEPENDENCIES Hist RIO Cint Core Gui Gpad)
      else()
        Set(DEPENDENCIES Hist RIO Cling Core Gui Gpad)
      endif()
      
      GENERATE_LIBRARY()

      Set(EXE_NAME runSolveDynEqRoot)
      Set(SRCS 
        run/runSolveDynEqRoot.cxx
      )
      #set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
      #set(LIBS_MAP_LAPACK "LAPACK")
      if (ROOT_FOUND_VERSION LESS 59999)
        Set(DEPENDENCIES breit_utils blas lapack gfortran ${LAPACK_LIBRARIES}
           Hist Graf Gpad RIO Cint Core handle_root_signal
           )
      else()
        Set(DEPENDENCIES breit_utils blas lapack gfortran ${LAPACK_LIBRARIES}
           Hist Graf Gpad RIO Cling Core handle_root_signal
           )
      endif()
      
      GENERATE_EXECUTABLE()



      Set(EXE_NAME runSolveDynEqRootRK)
      Set(SRCS 
        run/runSolveDynEqRoot_RK.cxx
      )
      #set(LIBS ${LIBS} ${LAPACK_LIBRARIES})
      #set(LIBS_MAP_LAPACK "LAPACK")
      if (ROOT_FOUND_VERSION LESS 59999)
        Set(DEPENDENCIES breit_utils blas lapack gfortran ${LAPACK_LIBRARIES}
           Hist Graf Gpad RIO Cint Core handle_root_signal
           )
      else()
        Set(DEPENDENCIES breit_utils blas lapack gfortran ${LAPACK_LIBRARIES}
           Hist Graf Gpad RIO Cling Core handle_root_signal
           )
      endif()
      GENERATE_EXECUTABLE()


  endif(ROOT_FOUND)

endif(LAPACK_FOUND)

