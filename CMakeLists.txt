cmake_minimum_required(VERSION 3.6)
project(freddi)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
        src/main.cpp
        src/arguments.hpp src/arguments.cpp
        src/gsl_const_cgsm.h
        src/nonlinear_diffusion.cpp
        src/nonlinear_diffusion.hpp
        src/opacity_related.cpp
        src/opacity_related.hpp
        src/orbit.cpp
        src/orbit.hpp
        src/spectrum.cpp
        src/spectrum.hpp src/constants.hpp
        src/freddi.cpp src/freddi.hpp
        src/output.cpp src/output.hpp
        src/unit_transfomation.hpp
)

SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "/usr/local/inlcude")
SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "/usr/local/lib")

FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
ADD_DEFINITIONS( "-DHAS_BOOST" )

add_executable(freddi ${SOURCE_FILES})

target_link_libraries(freddi LINK_PUBLIC ${Boost_LIBRARIES})
