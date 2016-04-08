cmake_minimum_required(VERSION 3.5)
project(GaussSeidelCPP_CL)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES GS_MainProcedure.cpp)
add_executable(GaussSeidelCPP_CL ${SOURCE_FILES})