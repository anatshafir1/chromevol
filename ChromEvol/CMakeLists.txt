# CMake script for chromevol
# Authors:
#   Anat Shafir (2022)
# Created: 06/12/2022

# Add all tests.
# Any .cpp file in test/ is considered to be a test.
# It will be compiled as a standalone program (must contain a main()).
# A test is considered to succeed if it returns EXIT_SUCCESS (usually 0).
# Tests are linked to the the shared library target.

file (GLOB cpp_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
foreach (cpp_file ${cpp_files})
  # Add each test (named as the filename without extension)
  get_filename_component (program_name ${cpp_file} NAME_WE)
  add_executable (${program_name} ${cpp_file})
  target_link_libraries (${program_name} ${PROJECT_NAME}-shared)
  set_target_properties (${program_name} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
  add_test (
    NAME ${program_name}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND ${program_name}
    )
endforeach (cpp_file)
