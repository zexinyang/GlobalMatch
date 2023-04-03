cmake_minimum_required(VERSION 3.12 FATAL_ERROR)

project(CSF LANGUAGES CXX)

set(package_name "CSF")

set(LIBRARY_OUTPUT_PATH ../lib)

# Zexin: to be compatible with PCL point clouds
find_package(PCL REQUIRED COMPONENTS common)
find_package(VTK REQUIRED NO_MODULE)

add_library(CSF STATIC
        c2cdist.cpp
        c2cdist.h
        Cloth.cpp
        Cloth.h
        Constraint.cpp
        Constraint.h
        CSF.cpp
        CSF.h
        Particle.cpp
        Particle.h
        point_cloud.cpp
        point_cloud.h
        Rasterization.cpp
        Rasterization.h
        Vec3.h
        XYZReader.cpp
        XYZReader.h
        )

target_include_directories(CSF PUBLIC ${PCL_INCLUDE_DIRS})
target_link_libraries(CSF PUBLIC ${PCL_LIBRARIES} ${VTK_LIBRARIES})

install(TARGETS CSF DESTINATION ../lib)