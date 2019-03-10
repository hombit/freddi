cmake_minimum_required(VERSION 3.0)
project(freddi)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -O3")

set(STATIC_LINKING FALSE CACHE BOOL "Build a static binary?")
if(STATIC_LINKING)
    set(CMAKE_EXE_LINKER_FLAGS "-static")
endif(STATIC_LINKING)

include_directories ("${PROJECT_SOURCE_DIR}/cpp/include")

find_package(Boost 1.57.0 COMPONENTS program_options unit_test_framework REQUIRED)
include_directories(${Boost_INCLUDE_DIR})
link_directories ( ${Boost_LIBRARY_DIRS} )

function(CREATE_EXE targ)
    set(TARGET freddi${targ})
    file(GLOB_RECURSE SOURCES cpp/main${targ}.cpp cpp/src/*.cpp cpp/include/*.hpp cpp/include/*.h)
    add_executable(${TARGET} ${SOURCES})
    target_compile_definitions(${TARGET} PUBLIC INSTALLPATHPREFIX="${CMAKE_INSTALL_PREFIX}")
    target_link_libraries(${TARGET} LINK_PUBLIC ${Boost_LIBRARIES})
    install(TARGETS ${TARGET} DESTINATION bin)
    install(FILES ${PROJECT_SOURCE_DIR}/freddi.ini DESTINATION etc)
endfunction()

CREATE_EXE("")
CREATE_EXE("-ns")

enable_testing()

add_executable(test_util cpp/test/util.cpp cpp/src/util.cpp)
if(NOT ${STATIC_LINKING})
    target_compile_definitions(test_util PUBLIC BOOST_TEST_DYN_LINK)
endif()
target_link_libraries(test_util LINK_PUBLIC ${Boost_LIBRARIES})
add_test(test_util test_util)
