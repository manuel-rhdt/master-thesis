cmake_minimum_required(VERSION 3.14)
project(Gillespie C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

include_directories(.)
include_directories(numtools)

add_compile_definitions(_BSD_SOURCE)

set(POW_LIBS "")
include(CheckLibraryExists)
check_library_exists(m pow "" LIBM)
if (LIBM)
    list(APPEND POW_LIBS "m")
endif ()

add_executable(Gillespie
        numtools/numtools.c
        numtools/vectools.c
        analyse.c
        main.c
        propagate.c
        ran.n.c
        read_comp_react.c
        stats.c propagate.h)

target_link_libraries(Gillespie PUBLIC ${POW_LIBS})