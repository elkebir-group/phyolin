cmake_minimum_required(VERSION 2.8)
project(phyolin)

set( LIBLEMON_ROOT "$ENV{HOME}/lemon" CACHE PATH "Additional search directory for lemon library" )

set( CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR} ${CMAKE_MODULE_PATH} )


ADD_DEFINITIONS(-DCPLEX -DIL_STD)
if( NOT CPLEX_INC_DIR OR NOT CPLEX_LIB_DIR OR NOT CONCERT_INC_DIR OR NOT CONCERT_LIB_DIR)
  set( CPLEX_AUTODETECTION_PREFIXES "$ENV{HOME}/ILOG;/Applications/;/opt/ibm/ILOG;$ENV{HOME}/Applications/IBM/ILOG" CACHE STRING "CPLEX automatic lookup directories" )
  foreach( CplexInstallPrefix ${CPLEX_AUTODETECTION_PREFIXES} )
    file( GLOB_RECURSE CplexLibraryPaths ${CplexInstallPrefix}/*/cplex/lib/*/libcplex.a )
    if( CplexLibraryPaths )
      list( LENGTH CplexLibraryPaths NumPath )
      if( NumPath GREATER 1 )
        message( WARNING "Found multiple CPLEX installations in ${CplexInstallPrefix}, automatically selecting the first one, please set {CPLEX|CONCERT}_{INC|LIB}_DIR manually if this is not the right selection" )
      else()
        message( "Found CPLEX's installation in ${CPlexInstallPrefix}" )
      endif()
      list( GET CplexLibraryPaths 0 CplexLibraryPath )

      set( MyRegex "${CplexInstallPrefix}/(.*)/cplex/?/lib/(.*)/libcplex.a" )
      string( REGEX MATCH "${MyRegex}" ValidatedCplexLibPath ${CplexLibraryPath} )
      if( NOT ValidatedCplexLibPath )
        message( WARNING "Unable to recognize ${CplexLibraryPath} as a valid CPLEX path" )
      else()
        string( REGEX REPLACE "${MyRegex}" "\\1" CplexHome ${ValidatedCplexLibPath} )
        string( REGEX REPLACE "${MyRegex}" "\\2" CplexArch ${ValidatedCplexLibPath} )

        set( CPLEX_INC_DIR   "${CplexInstallPrefix}/${CplexHome}/cplex/include"            CACHE PATH "CPLEX's include dir" )
        set( CP_INC_DIR   "${CplexInstallPrefix}/${CplexHome}/cpoptimizer/include"            CACHE PATH "CPOpt include dir" )
        set( CPLEX_LIB_DIR   "${CplexInstallPrefix}/${CplexHome}/cplex/lib/${CplexArch}"   CACHE PATH "CPLEX's library dir" )
        set( CP_LIB_DIR   "${CplexInstallPrefix}/${CplexHome}/cpoptimizer/lib/${CplexArch}"   CACHE PATH "CPopts library dir" )
        set( CONCERT_INC_DIR "${CplexInstallPrefix}/${CplexHome}/concert/include"          CACHE PATH "Concert's include dir" )
        set( CONCERT_LIB_DIR "${CplexInstallPrefix}/${CplexHome}/concert/lib/${CplexArch}" CACHE PATH "Concert's library dir" )
        
        break()
      endif()
    else()
      message( "Unable to find CPLEX and Concert installation in ${CplexInstallPrefix}" )
    endif()
  endforeach()

  if( NOT CPLEX_INC_DIR OR NOT CPLEX_LIB_DIR OR NOT CONCERT_INC_DIR OR NOT CONCERT_LIB_DIR )
    message( "Unable to autodeted CPLEX's and Concert's installation directories, please set {CPLEX|CONCERT}_{INC|LIB}_DIR manually" )
  else()
    message( "Autodetected CPLEX's include dir   : ${CPLEX_INC_DIR}" )
    message( "Autodetected CPLEX's library dir   : ${CPLEX_LIB_DIR}" )
    message( "Autodetected Concert's include dir : ${CONCERT_INC_DIR}" )
    message( "Autodetected Concert's library dir : ${CONCERT_LIB_DIR}" )
  endif()
endif()  

include_directories( ${CPLEX_INC_DIR} ${CONCERT_INC_DIR} "${LIBLEMON_ROOT}/include" ${CP_INC_DIR} "src" ${Boost_INCLUDE_DIRS} )
link_directories( ${CPLEX_LIB_DIR} ${CONCERT_LIB_DIR} ${CP_LIB_DIR} "${LIBLEMON_ROOT}/lib" )

set (phyolin_src
  src/phyolinCP.cpp
  src/singlecellmatrix.cpp
)

set (phyolin_hdr
  src/phyolinCP.h
  src/singlecellmatrix.h
)



set( Boost_USE_STATIC_LIBS ON )
#set( Boost_USE_SHARED_LIBS ON )
  find_package( Boost REQUIRED )
  MESSAGE( STATUS "found boost:"
  "\nINCLUDE: ${Boost_INCLUDE_DIRS}"
  "\nLIB: ${Boost_LIBRARIES}"
)
set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wno-unused-local-typedef" )




set( CommonLibs
  m
  pthread
  ${Boost_LIBRARIES}
  emon
  dl
)


add_executable( phyolin ${phyolin_src} ${phyolin_hdr} )
target_include_directories( phyolin PUBLIC "src" "${LIBLEMON_ROOT}/include" ${Boost_INCLUDE_DIRS} ${CPLEX_INC_DIR} ${CONCERT_INC_DIR} ${CP_INC_DIR})
target_link_libraries(phyolin cp  concert cp cplex  ${CommonLibs} )

