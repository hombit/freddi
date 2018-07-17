cmake_minimum_required(VERSION 3.6)
project(freddi)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

file(GLOB_RECURSE SOURCES cpp/*.cpp cpp/*.h cpp/*.hpp)

SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} "/usr/local/inlcude")
SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} "/usr/local/lib")

FIND_PACKAGE(Boost COMPONENTS program_options REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
ADD_DEFINITIONS( "-DHAS_BOOST" )

add_executable(freddi ${SOURCES})

target_link_libraries(freddi LINK_PUBLIC ${Boost_LIBRARIES})
