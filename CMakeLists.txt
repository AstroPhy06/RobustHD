cmake_minimum_required(VERSION 3.16)

# set the project name
project(RobustHD)

set(CMAKE_CXX_STANDARD 17)

find_package (OpenCV REQUIRED )
include_directories( ${OpenCV_INCLUDE_DIRS} )

find_package(OpenMP)

# add the executable
add_executable(RobustHD src/main.cpp src/grid.hpp src/solveur.hpp src/solveur.cpp src/physics.hpp src/scheme.hpp src/scheme.cpp)
target_link_libraries(RobustHD ${OpenCV_LIBS} "${OpenMP_CXX_FLAGS}")
target_compile_options(RobustHD PRIVATE "${OpenMP_CXX_FLAGS}")
