cmake_minimum_required(VERSION 3.14)
project(Gillespie CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories(.)
include_directories(numtools)
add_compile_definitions(_BSD_SOURCE)

add_subdirectory(numtools)

set(JSON_BuildTests OFF CACHE INTERNAL "")
set(JSON_Install OFF CACHE INTERNAL "")
add_subdirectory(json)

add_executable(Gillespie
        main.cpp Simulation.cpp Simulation.hh Run.cpp Run.hh Block.cpp Block.hh Trajectory.cpp Trajectory.hh Reaction.cpp Reaction.hh)

target_link_libraries(Gillespie PRIVATE numtools)
target_link_libraries(
        Gillespie PRIVATE nlohmann_json::nlohmann_json)
