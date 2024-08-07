project(DREndcapTubes)

file(GLOB headers
  ${PROJECT_SOURCE_DIR}/include/*.h
)

dd4hep_add_plugin(DREndcapTubes SOURCES src/*.cpp USES
  DD4hep::DDCore
  DD4hep::DDCond
  ROOT::Core
  ROOT::Geom
  ROOT::GenVector
  ROOT::MathCore
)

target_include_directories(DREndcapTubes PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
)

# Link the Utils library
#target_link_libraries(DREndcapTubes PRIVATE DRutils)

set_target_properties(DREndcapTubes PROPERTIES PUBLIC_HEADER "${headers}")

install(TARGETS DREndcapTubes
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}" COMPONENT dev
)

install(DIRECTORY compact DESTINATION ${CMAKE_INSTALL_DATADIR})

dd4hep_configure_scripts( DREndcapTubes DEFAULT_SETUP WITH_TESTS )
