# bin/CMakeLists.txt

set(BIN_SOURCES
  test_deriv_costg.cpp
  eom.cpp
  eomboost.cpp
)

# Process each source file and create an executable
foreach(SOURCE_FILE IN LISTS BIN_SOURCES)
  # Get the filename without an extension
  get_filename_component(EXECUTABLE_NAME ${SOURCE_FILE} NAME_WE)

  # Define the executable
  add_executable(${EXECUTABLE_NAME} ${SOURCE_FILE})

  # Include the directory for headers
  # target_include_directories(${EXECUTABLE_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

  # Link the executable to the required libraries
  target_link_libraries(${EXECUTABLE_NAME} PRIVATE 
    integrator 
    yaml-cpp
    dso::sp3
    ${PROJECT_DEPENDENCIES}
    ${CSPICE_LIB}
  )
endforeach()
