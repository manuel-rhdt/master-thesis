cmake_minimum_required(VERSION 3.14)
project(numtools)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

add_library(numtools STATIC
        numtools.c
        numtools.h
        vectools.c
        vectools.h)