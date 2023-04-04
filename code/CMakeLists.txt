cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(GlobalMatch)

# Set C++ standard
set(CMAKE_CXX_STANDARD 14)

# Find packages
find_package(VTK REQUIRED NO_MODULE)

find_package(PCL REQUIRED COMPONENTS
        common io features filters segmentation
        sample_consensus kdtree registration
        CONFIG
        )

find_package(OpenMP)
if (OpenMP_FOUND)
    # Add OpenMP flags to all targets
    add_compile_options(${OpenMP_CXX_FLAGS})
    link_libraries(${OpenMP_CXX_FLAGS})
endif ()

# Add subdirectory for the CSF library
add_subdirectory(3rd_party/csf)

# Set project sources
set(PROJECT_SOURCES
        3rd_party/voxel_grid/voxel_grid_fix.cpp
        3rd_party/octree/octree_extract_clusters.cpp
        global_match/main.cpp
        global_match/stem_mapping.cpp
        global_match/stem_matching.cpp
        )

# Add executable target
add_executable(GlobalMatch ${PROJECT_SOURCES})

# Add include directories to target
target_include_directories(GlobalMatch PUBLIC
        ${PCL_INCLUDE_DIRS} # PCL include directories
        3rd_party/cc
        3rd_party/csf
        3rd_party/octree
        3rd_party/voxel_grid
        )

# Link libraries to target
target_link_libraries(GlobalMatch PUBLIC
        CSF # CSF library
        ${PCL_LIBRARIES} # PCL libraries
        ${VTK_LIBRARIES} # VTK libraries
        )