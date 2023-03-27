project(CSF LANGUAGES CXX)

set(package_name "CSF")

set(LIBRARY_OUTPUT_PATH ../lib)

# Zexin: to be compatible with PCL point clouds
find_package(PCL REQUIRED COMPONENTS common)
find_package(VTK COMPONENTS CommonCore CommonDataModel)

add_library(CSF STATIC
        ${CMAKE_CURRENT_LIST_DIR}/c2cdist.cpp
        ${CMAKE_CURRENT_LIST_DIR}/c2cdist.h
        ${CMAKE_CURRENT_LIST_DIR}/Cloth.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Cloth.h
        ${CMAKE_CURRENT_LIST_DIR}/Constraint.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Constraint.h
        ${CMAKE_CURRENT_LIST_DIR}/CSF.cpp
        ${CMAKE_CURRENT_LIST_DIR}/CSF.h
        ${CMAKE_CURRENT_LIST_DIR}/Particle.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Particle.h
        ${CMAKE_CURRENT_LIST_DIR}/point_cloud.cpp
        ${CMAKE_CURRENT_LIST_DIR}/point_cloud.h
        ${CMAKE_CURRENT_LIST_DIR}/Rasterization.cpp
        ${CMAKE_CURRENT_LIST_DIR}/Rasterization.h
        ${CMAKE_CURRENT_LIST_DIR}/Vec3.h
        ${CMAKE_CURRENT_LIST_DIR}/XYZReader.cpp
        ${CMAKE_CURRENT_LIST_DIR}/XYZReader.h
        )

include_directories(${PCL_INCLUDE_DIRS})
link_directories(${PCL_LIBRARY_DIRS})
add_definitions(${PCL_DEFINITIONS})

target_link_libraries(CSF PUBLIC ${VTK_LIBRARIES})

install(TARGETS CSF DESTINATION ../lib)