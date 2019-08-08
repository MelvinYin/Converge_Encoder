cmake_minimum_required(VERSION 3.14)
project(conv_fasta_encoder)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_INCLUDE_CURRENT_DIR 1)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(conv_fasta_encoder main.cpp)