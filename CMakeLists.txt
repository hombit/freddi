cmake_minimum_required(VERSION 3.12)
project(freddi)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

include_directories("${PROJECT_SOURCE_DIR}/cpp/include")

function(CREATE_EXE targ)
    set(TARGET freddi${targ}-exe)

    set(STATIC_LINKING FALSE CACHE BOOL "Build a static binary?")
    if(STATIC_LINKING)
        set(CMAKE_EXE_LINKER_FLAGS "-static")
    endif(STATIC_LINKING)


    file(GLOB_RECURSE SOURCES cpp/main${targ}.cpp cpp/src/*.cpp cpp/include/*.hpp cpp/include/*.h)
    add_executable(${TARGET} ${SOURCES})

    find_package(Boost 1.57.0 COMPONENTS program_options REQUIRED)
    target_include_directories(${TARGET} PUBLIC ${Boost_INCLUDE_DIR})

    target_compile_definitions(${TARGET} PUBLIC INSTALLPATHPREFIX="${CMAKE_INSTALL_PREFIX}")
    target_link_libraries(${TARGET} LINK_PUBLIC ${Boost_LIBRARIES})
    install(TARGETS ${TARGET} DESTINATION bin)
    install(FILES ${PROJECT_SOURCE_DIR}/freddi.ini DESTINATION etc)
endfunction()

CREATE_EXE("")
CREATE_EXE("-ns")

set(TARGET _freddi)
find_package(Python3 COMPONENTS Development Interpreter REQUIRED)
set(PYVERSION ${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR})
find_package(Boost 1.64.0 COMPONENTS program_options python${PYVERSION} numpy${PYVERSION} REQUIRED)
file(GLOB_RECURSE SOURCES cpp/src/*.cpp cpp/include/*.hpp cpp/include/*.h cpp/pywrap/*.cpp cpp/pywrap/*.hpp)
Python3_add_library(${TARGET} MODULE ${SOURCES})
if (MODULE_OUTPUT_NAME)
    set_property(TARGET ${TARGET} PROPERTY LIBRARY_OUTPUT_NAME ${MODULE_OUTPUT_NAME})
endif()
if (MODULE_OUTPUT_SUFFIX)
    set_property(TARGET ${TARGET} PROPERTY SUFFIX ${MODULE_OUTPUT_SUFFIX})
endif()
execute_process(COMMAND ${Python3_EXECUTABLE} -c "import numpy; print(numpy.get_include(), end='')" OUTPUT_VARIABLE NumPy_INCLUDE_DIR)
target_include_directories(${TARGET} PUBLIC ${Python3_INCLUDE_DIRS} ${NumPy_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
target_link_libraries(${TARGET} LINK_PUBLIC ${Boost_LIBRARIES} ${Python3_LIBRARIES})
