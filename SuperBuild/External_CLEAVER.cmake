
set(proj CLEAVER)

# Set dependency list
set(${proj}_DEPENDENCIES "")

# Include dependent projects if any
ExternalProject_Include_Dependencies(${proj} PROJECT_VAR proj DEPENDS_VAR ${proj}_DEPENDENCIES)

# Sanity checks
if(DEFINED ${proj}_INCLUDE_DIR AND NOT EXISTS ${${proj}_INCLUDE_DIR})
  message(FATAL_ERROR "${proj}_INCLUDE_DIR variable is defined but corresponds to nonexistent directory")
endif()
if(DEFINED ${proj}_LIBRARY AND NOT EXISTS ${${proj}_LIBRARY})
  message(FATAL_ERROR "${proj}_LIBRARY variable is defined but corresponds to nonexistent file")
endif()

if(NOT DEFINED ${proj}_INCLUDE_DIR OR NOT DEFINED ${proj}_LIBRARY)

  set(EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS )

  if(NOT CMAKE_CONFIGURATION_TYPES)
    list(APPEND EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS
      -DCMAKE_BUILD_TYPE:STRING=${CMAKE_BUILD_TYPE}
      )
  endif()

  set(EP_SOURCE_DIR ${CMAKE_BINARY_DIR}/${proj})
  set(EP_BINARY_DIR ${CMAKE_BINARY_DIR}/${proj}-build)

  ExternalProject_add(${proj}
    GIT_REPOSITORY https://github.com/SCIInstitute/Cleaver2.git
    GIT_TAG bb1d7e0604171b52f4e2f3b64178860dde0e4a58
    SOURCE_DIR ${EP_SOURCE_DIR}
    BINARY_DIR ${EP_BINARY_DIR}
    CMAKE_CACHE_ARGS
      # Options
      -DBUILD_CLI:BOOL=OFF
      -DBUILD_EXTRAS:BOOL=OFF
      -DBUILD_GUI:BOOL=OFF
      -DBUILD_SHARED_LIBS:BOOL=OFF
      -DBUILD_TESTING:BOOL=OFF
      -DTeem_BZIP2:BOOL=OFF
      -DUSE_GCOV:BOOL=OFF
      -DUSE_ITK:BOOL=OFF
      -DCOPY_HEADER_FILES:BOOL=ON
      # Install directories
      # NA
      ${EXTERNAL_PROJECT_OPTIONAL_CMAKE_CACHE_ARGS}
    INSTALL_COMMAND ""
    USES_TERMINAL_DOWNLOAD 1
    USES_TERMINAL_UPDATE 1
    USES_TERMINAL_CONFIGURE 1
    USES_TERMINAL_BUILD 1
    USES_TERMINAL_INSTALL 1
    )

  set(${proj}_INCLUDE_DIR ${EP_BINARY_DIR}/include)

  if(WIN32)
    set(${proj}_LIBRARY ${EP_BINARY_DIR}/lib/libcleaver.lib)
  else()
    set(${proj}_LIBRARY ${EP_BINARY_DIR}/lib/libcleaver.a)
  endif()

endif()

mark_as_superbuild(
  VARS
    ${proj}_INCLUDE_DIR:PATH
    ${proj}_LIBRARY:FILEPATH
  )

ExternalProject_Message(${proj} "${proj}_INCLUDE_DIR:${${proj}_INCLUDE_DIR}")
ExternalProject_Message(${proj} "${proj}_LIBRARY:${${proj}_LIBRARY}")
