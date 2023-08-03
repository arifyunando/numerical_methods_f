cmake_minimum_required(VERSION 3.15)

# set the project name and version
project(... VERSION 1.0)
enable_language(fortran)

# add the libraries

# add the executable
add_executable(${PROJECT_NAME} ....f90)

# add the install targets
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/bin)

# enable testing