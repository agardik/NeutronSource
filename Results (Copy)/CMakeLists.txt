cmake_minimum_required(VERSION 3.16)
project(TreeReader)

# Find ROOT
find_package(ROOT REQUIRED)

# Include ROOT libraries
include(${ROOT_USE_FILE})

# Add executables
add_executable(DoseCalculation DoseCalculation.cxx TreeReader.cxx)
target_link_libraries(DoseCalculation ${ROOT_LIBRARIES})
