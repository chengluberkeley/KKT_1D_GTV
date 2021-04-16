cmake_minimum_required (VERSION 3.0)
project(KKT)

if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()
set(CMAKE_CXX_FLAGS "-std=c++11 -O3")

include_directories(
  ${PROJECT_SOURCE_DIR}
  ${PROJECT_SOURCE_DIR}/KKT
  ${PROJECT_SOURCE_DIR}/Data
  ${PROJECT_SOURCE_DIR}/Profiles
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)

file(GLOB lib_srcs
  ${PROJECT_SOURCE_DIR}/Profiles/*.cpp
  ${PROJECT_SOURCE_DIR}/Data/*.cpp
  ${PROJECT_SOURCE_DIR}/KKT/*.cpp
)

add_executable(kkt_main ${PROJECT_SOURCE_DIR}/main.cpp
  ${lib_srcs}
)

