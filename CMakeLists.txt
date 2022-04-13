# Cleaver library
find_library(CLEAVER_LIBRARY NAMES cleaver)
if(NOT CLEAVER_LIBRARY)
  message(FATAL_ERROR "Could not find Cleaver library")
endif()
message(STATUS "Found Cleaver library at: ${CLEAVER_LIBRARY}")

# Cleaver include directory
find_path(CLEAVER_INCLUDE_DIR cleaver/Cleaver.h)
if(NOT CLEAVER_INCLUDE_DIR)
  message(FATAL_ERROR "Could not find Cleaver include directory")
endif()
message(STATUS "Found Cleaver include directory: ${CLEAVER_INCLUDE_DIR}")

# JsonCpp library
if(NOT DEFINED CLEAVER_JsonCpp_LIBRARY)
  set(CLEAVER_JsonCpp_LIBRARY ${CLEAVER_LIBRARY})
  cmake_path(REMOVE_FILENAME CLEAVER_JsonCpp_LIBRARY)
  if(WIN32)
    cmake_path(APPEND CLEAVER_JsonCpp_LIBRARY "jsoncpp" "jsoncpp.lib")
  elseif(APPLE)
    # TODO
  else()
    cmake_path(APPEND CLEAVER_JsonCpp_LIBRARY "libjsoncpp.a")
  endif()
endif()
if(NOT CLEAVER_JsonCpp_LIBRARY)
  message(FATAL_ERROR "Could not find Cleaver JsonCpp library")
endif()
message(STATUS "Found Cleaver JsonCpp library at: ${CLEAVER_JsonCpp_LIBRARY}")

# JsonCpp include directory
if(NOT DEFINED CLEAVER_JsonCpp_INCLUDE_DIR)
  set(CLEAVER_JsonCpp_INCLUDE_DIR ${CLEAVER_INCLUDE_DIR})
endif()
if(NOT CLEAVER_JsonCpp_INCLUDE_DIR)
  message(FATAL_ERROR "Could not find Cleaver JsonCpp include directory")
endif()
message(STATUS "Found Cleaver JsonCpp include directory at: ${CLEAVER_JsonCpp_INCLUDE_DIR}")

# Cleaver target
add_library(Cleaver::Cleaver UNKNOWN IMPORTED)
set_target_properties(Cleaver::Cleaver PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${CLEAVER_INCLUDE_DIR}"
  )
set_property(TARGET Cleaver::Cleaver APPEND PROPERTY
  IMPORTED_LOCATION "${CLEAVER_LIBRARY}"
  )

# JsonCpp target
add_library(Cleaver::JsonCpp UNKNOWN IMPORTED)
set_target_properties(Cleaver::JsonCpp PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${CLEAVER_JsonCpp_INCLUDE_DIR}"
  )
set_property(TARGET Cleaver::JsonCpp APPEND PROPERTY
  IMPORTED_LOCATION "${CLEAVER_JsonCpp_LIBRARY}"
  )

# VTK module
set(classes
  vtkCleaverImageToUnstructuredGridFilter
  )

vtk_module_add_module(VTK::Cleaver
  CLASSES ${classes}
  )
vtk_module_link(VTK::Cleaver
  PUBLIC
    Cleaver::Cleaver
    Cleaver::JsonCpp
  )
