set(PROJECT_NAME vio_init_sys)

PROJECT(${PROJECT_NAME})
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

set(CMAKE_BUILD_TYPE Release) # Release, RelWithDebInfo , Debug

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fext-numeric-literals")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS}  -Wall  -O3 -march=native ")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall  -O3 -march=native")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake_modules/")

add_definitions("-D SOPHUS_USE_BASIC_LOGGING")
add_definitions(-DBOOST_ERROR_CODE_HEADER_ONLY)

find_package(Eigen3 REQUIRED)
find_package(OpenCV REQUIRED)
find_package(Pangolin REQUIRED)
find_package(Ceres 1.14.0 REQUIRED)
find_package(Boost REQUIRED)
include_directories(
        ${EIGEN3_INCLUDE_DIR}
        ${OPENCV_INCLUDE_DIR}
        ${CERES_INCLUDE_DIRS}
        ${Pangolin_INCLUDE_DIRS}
        ${Boost_INCLUDE_DIR}
        ${PROJECT_SOURCE_DIR}/thirdparty/sophus
        ${PROJECT_SOURCE_DIR}/thirdparty/spectra
        ${PROJECT_SOURCE_DIR}/src
        "usr/local/inclue/eigen3"
)

add_library(${PROJECT_NAME} SHARED
        #  simulation code
        src/param.h
        src/param.cpp
        src/utilities.h
        src/utilities.cpp
        src/imu.h
        src/imu.cpp
        src/viewer.h
        src/viewer.cpp
        src/linearCoplanarError.cpp
        src/random.cpp
        src/evaluation.cc
        src/coplanar/point_projection_factor.cc
        src/coplanar/pose_local_parameterization.cc
        )


target_link_libraries(${PROJECT_NAME}
        ${Pangolin_LIBRARIES}
        ${OpenCV_LIBS}
        ${CERES_LIBRARIES}
        ${Boost_LIBRARIES}
        ${Sophus_LIBRARIES}
        )

ADD_EXECUTABLE(data_gen main/gener_alldata.cpp)
TARGET_LINK_LIBRARIES(data_gen ${PROJECT_NAME})