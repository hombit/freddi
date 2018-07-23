cmake_minimum_required(VERSION 3.0)
project(freddi)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

set(STATIC_LINKING FALSE CACHE BOOL "Build a static binary?")
if(STATIC_LINKING)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
endif(STATIC_LINKING)

find_package(Boost 1.57.0 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
if(Boost_PROGRAM_OPTIONS_FOUND)
    add_definitions("-DHAS_BOOST_PROGRAM_OPTIONS")
endif(Boost_PROGRAM_OPTIONS_FOUND)

set(TARGET freddi)
file(GLOB_RECURSE SOURCES cpp/main.cpp cpp/src/*.cpp cpp/include/*.hpp cpp/include/*.h)
include_directories ("${PROJECT_SOURCE_DIR}/cpp/include")
add_executable(${TARGET} ${SOURCES})
target_compile_definitions(${TARGET} PUBLIC INSTALLPATHPREFIX="${CMAKE_INSTALL_PREFIX}")
target_link_libraries(${TARGET} LINK_PUBLIC ${Boost_LIBRARIES})
install(TARGETS ${TARGET} DESTINATION bin)
install(FILES ${PROJECT_SOURCE_DIR}/freddi.ini DESTINATION etc)

