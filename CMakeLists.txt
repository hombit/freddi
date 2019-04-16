cmake_minimum_required(VERSION 3.1)
message("CMake version: ${CMAKE_VERSION}")
project(freddi)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

include_directories("${PROJECT_SOURCE_DIR}/cpp/include")

find_package(Boost 1.57.0 COMPONENTS program_options REQUIRED)

function(CREATE_EXE targ)
    set(TARGET freddi${targ}-exe)
    set(EXE freddi${targ})

    set(STATIC_LINKING FALSE CACHE BOOL "Build a static binary?")
    if(STATIC_LINKING)
        set(Boost_USE_STATIC_LIBS "ON")
    endif(STATIC_LINKING)


    file(GLOB_RECURSE SOURCES cpp/main${targ}.cpp cpp/src/*.cpp cpp/include/*.hpp cpp/include/*.h)
    add_executable(${TARGET} ${SOURCES})
    set_property(TARGET ${TARGET} PROPERTY OUTPUT_NAME ${EXE})

    target_include_directories(${TARGET} PUBLIC ${Boost_INCLUDE_DIR})

    target_compile_definitions(${TARGET} PUBLIC INSTALLPATHPREFIX="${CMAKE_INSTALL_PREFIX}")
    target_link_libraries(${TARGET} ${Boost_LIBRARIES})
    install(TARGETS ${TARGET} DESTINATION bin)
    install(FILES ${PROJECT_SOURCE_DIR}/freddi.ini DESTINATION etc)
endfunction()

CREATE_EXE("")
CREATE_EXE("-ns")


function(CREATE_PY_MODULE)
    find_package(Python3 COMPONENTS Development Interpreter)
    foreach(Boost_Python3_VERSION_TAG
            "${Python3_VERSION_MAJOR}"
            "${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}"
            "-py${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}"
            "-${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}"
            "${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}"
            "-${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}")
        find_package(Boost QUIET COMPONENTS "python${Boost_Python3_VERSION_TAG}")
        if(${Boost_FOUND})
            message("Boost_Python3_VERSION_TAG: ${Boost_Python3_VERSION_TAG}")
            find_package(Boost COMPONENTS "python${Boost_Python3_VERSION_TAG}" "numpy${Boost_Python3_VERSION_TAG}")
            break()
        endif()
    endforeach()
    if (${Python3_FOUND} AND ${Boost_FOUND})
        set(TARGET _freddi)
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
    endif()
endfunction(CREATE_PY_MODULE)

if (NOT NO_PYTHON_MODULE)
    CREATE_PY_MODULE()
endif()
