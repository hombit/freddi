cmake_minimum_required(VERSION 3.6)
project(freddi)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -DINSTALLPATHPREFIX=''")

set(SOURCE_FILES
        main.cpp
    arguments.hpp arguments.cpp
    gsl_const_cgsm.h
    nonlinear_diffusion.cpp
    nonlinear_diffusion.hpp
    opacity_related.cpp
    opacity_related.hpp
    orbit.cpp
    orbit.hpp
    spectrum.cpp
    spectrum.hpp constants.hpp
    freddi.cpp freddi.hpp
    output.cpp output.hpp
    unit_transfomation.hpp
)

SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "/usr/local/inlcude")
SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "/usr/local/lib")

FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
ADD_DEFINITIONS( "-DHAS_BOOST" )

add_executable(freddi ${SOURCE_FILES})

target_link_libraries(freddi LINK_PUBLIC ${Boost_LIBRARIES})