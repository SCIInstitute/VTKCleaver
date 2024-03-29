include(FetchContent)

set(proj VTKExternalModule)
if (FETCH_${proj}_INSTALL_LOCATION)
  # The install location can be specified
  set(EP_SOURCE_DIR "${FETCH_${proj}_INSTALL_LOCATION}")
else()
  set(EP_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj})
endif()

FetchContent_Populate(${proj}
  SOURCE_DIR     ${EP_SOURCE_DIR}
  GIT_REPOSITORY https://github.com/KitwareMedical/VTKExternalModule.git
  GIT_TAG        c1906cf121e34b6391a91c2fffc448eca402a6cc
  QUIET
  )

message(STATUS "Remote - ${proj} [OK]")

set(VTKExternalModule_SOURCE_DIR ${EP_SOURCE_DIR})
message(STATUS "Remote - VTKExternalModule_SOURCE_DIR:${VTKExternalModule_SOURCE_DIR}")
