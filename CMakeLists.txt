cmake_minimum_required(VERSION 3.8)
project(freehand3d_US_reconstruction)

enable_language(CUDA)

#enable c++11
set(CMAKE_CXX_STANDARD 14)

#fine CUDA:
find_package(CUDA)

#set YAML:
find_package(yaml-cpp REQUIRED)

include_directories(freehand3d_US_reconstruction ${CUDA_INCLUDE_DIRS} ${PROJECT_SOURCE_DIR})

add_executable(
    FreehandReconstruction 
    utilities/ultrasound_frames.cpp
    utilities/volumetric.cpp
    utilities/utilities.cpp
    src/FreehandReconstruction.cpp
    src/freehand_reconstruction_GPU.cu
    src/main.cpp
)

target_link_libraries(
    FreehandReconstruction 
    ${CUDA_LIBRARIES}
    yaml-cpp
)
