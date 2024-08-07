# -------------------------------------------------------------------------
# File: CMakeLists.txt
#
# Purpose: CMake file to compile the example MyDetector
#
# Author: Lorenzo Pezzotti, wildly using dd4hep examples
#
# Created: 2/5/2024
# -------------------------------------------------------------------------

#==========================================================================
#  AIDA Detector description implementation 
#--------------------------------------------------------------------------
# Copyright (C) Organisation europeenne pour la Recherche nucleaire (CERN)
# All rights reserved.
#
# For the licensing terms see $DD4hepINSTALL/LICENSE.
# For the list of contributors see $DD4hepINSTALL/doc/CREDITS.
#
#==========================================================================

cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

IF(NOT TARGET DD4hep::DDCore)
  find_package ( DD4hep REQUIRED )
  include ( ${DD4hep_DIR}/cmake/DD4hep.cmake )
  include ( ${DD4hep_DIR}/cmake/DD4hepBuild.cmake )
  dd4hep_configure_output()
ENDIF()

dd4hep_set_compiler_flags()
#==========================================================================
if(TARGET XercesC::XercesC)
  SET(OPT_XERCESC XercesC::XercesC)
endif()

#Add plugin code with geometry for MyDetector
dd4hep_add_plugin(DREndcapTubes SOURCES ../src/*.cpp
    USES DD4hep::DDCore
         DD4hep::DDRec
         DD4hep::DDCond
         ROOT::Core
         ROOT::Geom
         ROOT::GenVector
         ${OPT_XERCESC}
  )

# Specify the include directory for the MyDetector target
target_include_directories(DREndcapTubes PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../include)

#Use this only if you do not want to use example from build directory
#install(TARGETS MyDetector LIBRARY DESTINATION lib)
#install(DIRECTORY compact DESTINATION examples/ClientTests)
#-----------------------------------------------------------------------------------
dd4hep_configure_scripts(DREndcapTubes DEFAULT_SETUP WITH_TESTS)
