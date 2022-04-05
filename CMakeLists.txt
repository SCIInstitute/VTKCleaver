# Cleaver external lib dependency
find_library(CLEAVER_LIBRARY NAMES cleaver)
if(NOT CLEAVER_LIBRARY)
  message(FATAL_ERROR "Could not find Cleaver library")
endif()
message(STATUS "Found Cleaver library at: ${CLEAVER_LIBRARY}")

# add additional libs
set(json_lib ${CLEAVER_LIBRARY})
cmake_path(REMOVE_FILENAME json_lib)
# message("path is \"${json_path}\"")
# cmake_path(GET json_path PARENT_PATH json_lib)
# cmake_path(REMOVE_FILENAME json_lib)
# message("Parent path is \"${json_lib}\"")
cmake_path(APPEND json_lib "jsoncpp" "jsoncpp.lib")

find_path(CLEAVER_INCLUDE_DIR cleaver/Cleaver.h)
if(NOT CLEAVER_INCLUDE_DIR)
  message(FATAL_ERROR "Could not find Cleaver include directory")
endif()
message(STATUS "Found Cleaver include directory: ${CLEAVER_INCLUDE_DIR}")

add_library(Cleaver::Cleaver UNKNOWN IMPORTED)
set_target_properties(Cleaver::Cleaver PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${CLEAVER_INCLUDE_DIR}"
  )
set_property(TARGET Cleaver::Cleaver APPEND PROPERTY
  IMPORTED_LOCATION "${CLEAVER_LIBRARY}"
  )
add_library(Cleaver::JSONCPP UNKNOWN IMPORTED)
set_property(TARGET Cleaver::JSONCPP APPEND PROPERTY
  IMPORTED_LOCATION "${json_lib}"
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
    Cleaver::JSONCPP
  )
