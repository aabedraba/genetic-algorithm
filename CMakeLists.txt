cmake_minimum_required(VERSION 2.8)
project(2-genetic-algorithms)

set(CMAKE_CXX_STANDARD 14)

add_executable(2-genetic-algorithms
        src/main.cpp
        src/Airport.cpp src/Airport.h
        src/Genetic.h src/Genetic.cpp
        src/Utils.h src/Utils.cpp)