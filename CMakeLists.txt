cmake_minimum_required(VERSION 3.14)
project(Gillespie C)

set(CMAKE_C_STANDARD 11)

include_directories(.)
include_directories(numtools)

set(POW_LIBS "")
include(CheckLibraryExists)
check_library_exists(m pow "" LIBM)
if (LIBM)
    list(APPEND POW_LIBS "m")
endif ()

add_executable(Gillespie
        numtools/numtools.c
        numtools/numtools.h
        numtools/vectools.c
        numtools/vectools.h
        analyse.c
        Gillespie.h
        main.c
        propagate.c
        ran.n.c
        read_comp_react.c
        stats.c
        stats.h)

target_link_libraries(Gillespie PUBLIC ${POW_LIBS})