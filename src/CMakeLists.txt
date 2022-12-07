# CMake script for ChromEvol lib
# Author:
#   Anat Shafir
# Created: 06/12/2022

# File list
set (CPP_FILES
  ChromosomeSubstitutionModel.cpp
  ChromEvolOptions.cpp
  UndirectedGraph.cpp
  ChromosomeNumberOptimizer.cpp
  ComputeChromosomeTransitionsExp.cpp
  ChromosomeNumberMng.cpp
  )

IF(BUILD_STATIC)
  # Build the static lib
  add_library (${PROJECT_NAME}-static STATIC ${CPP_FILES})
  target_include_directories (${PROJECT_NAME}-static PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
    )
  set_target_properties (${PROJECT_NAME}-static PROPERTIES OUTPUT_NAME ${PROJECT_NAME})
  target_link_libraries (${PROJECT_NAME}-static ${BPP_LIBS_STATIC} Eigen3::Eigen)
ENDIF()

# Build the shared lib
add_library (${PROJECT_NAME}-shared SHARED ${CPP_FILES})
target_include_directories (${PROJECT_NAME}-shared PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )
set_target_properties (${PROJECT_NAME}-shared
  PROPERTIES OUTPUT_NAME ${PROJECT_NAME}
  MACOSX_RPATH 1
  VERSION ${${PROJECT_NAME}_VERSION}
  SOVERSION ${${PROJECT_NAME}_VERSION_MAJOR}
  )
target_link_libraries (${PROJECT_NAME}-shared ${BPP_LIBS_SHARED} Eigen3::Eigen)

# Install libs and headers
IF(BUILD_STATIC)
install (
  TARGETS ${PROJECT_NAME}-static
  EXPORT ${PROJECT_NAME}-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
ENDIF()

install (
  TARGETS ${PROJECT_NAME}-shared
  EXPORT ${PROJECT_NAME}-targets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_LIBDIR}
  )
install (
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Bpp
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.h"
  )
# Install cmake file representing targets
install (EXPORT ${PROJECT_NAME}-targets DESTINATION ${cmake-package-location})