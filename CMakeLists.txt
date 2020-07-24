# Dependencies: libjpeg-devel, libpng-devel
cmake_minimum_required(VERSION 3.10)
project(RayTrace)
find_package(OpenMP REQUIRED)
find_package(JPEG REQUIRED)
find_package(PNG REQUIRED)
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fopenmp")

add_definitions("-DGLM_FORCE_RADIANS")

add_library(imageIO STATIC
        lib/imageIO/imageIO.cpp)

add_executable(${PROJECT_NAME}
        src/raytrace.cpp
        src/raytrace.h
)

target_link_libraries(${PROJECT_NAME} PRIVATE
        imageIO
        jpeg
        png
        OpenMP::OpenMP_CXX
)
target_include_directories(${PROJECT_NAME} PRIVATE
        lib/imageIO
)