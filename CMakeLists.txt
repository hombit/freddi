cmake_minimum_required(VERSION 3.0)
project(freddi)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

set(STATIC_LINKING FALSE CACHE BOOL "Build a static binary?")
if(STATIC_LINKING)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
endif(STATIC_LINKING)

file(GLOB_RECURSE SOURCES cpp/*.cpp cpp/*.h cpp/*.hpp)

find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

set(TARGET freddi)
add_executable(${TARGET} ${SOURCES})
target_compile_definitions(${TARGET} PUBLIC INSTALLPATHPREFIX="${CMAKE_INSTALL_PREFIX}")
target_link_libraries(${TARGET} LINK_PUBLIC ${Boost_LIBRARIES})
install(TARGETS ${TARGET} DESTINATION bin)
install(FILES ${PROJECT_SOURCE_DIR}/freddi.ini DESTINATION etc)

