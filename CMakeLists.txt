cmake_minimum_required(VERSION 3.1)
message("CMake version: ${CMAKE_VERSION}")
project(freddi)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

include_directories("${PROJECT_SOURCE_DIR}/cpp/include")


function(CREATE_SHARED)
    set(TARGET freddi-dll)
    find_package(Boost 1.57.0 COMPONENTS program_options filesystem REQUIRED)

    file(GLOB_RECURSE SOURCES
            cpp/src/*.cpp cpp/include/*.hpp cpp/include/*.h
            cpp/src/ns/*.cpp cpp/include/ns/*.hpp)
    add_library(${TARGET} SHARED ${SOURCES})
    set_property(TARGET ${TARGET} PROPERTY OUTPUT_NAME freddi)

    target_include_directories(${TARGET} PUBLIC ${Boost_INCLUDE_DIR})
    target_compile_definitions(${TARGET} PUBLIC INSTALLPATHPREFIX="${CMAKE_INSTALL_PREFIX}")
    target_link_libraries(${TARGET} ${Boost_LIBRARIES})
    install(TARGETS ${TARGET} DESTINATION lib)
endfunction()

CREATE_SHARED()


function(CREATE_EXE targ)
    set(TARGET freddi${targ}-exe)
    set(EXE freddi${targ})

    set(STATIC_LINKING FALSE CACHE BOOL "Build a static binary?")
    if(STATIC_LINKING)
        set(Boost_USE_STATIC_LIBS "ON")
    endif(STATIC_LINKING)

    find_package(Boost 1.57.0 COMPONENTS program_options filesystem REQUIRED)

    file(GLOB SOURCES
            cpp/main${targ}.cpp
            cpp/src/*.cpp cpp/include/*.hpp cpp/include/*.h)
    set(SOURCES_NS "")
    if (targ STREQUAL "-ns")
        file(GLOB SOURCES_NS cpp/src/ns/*.cpp cpp/include/ns/*.hpp)
    endif()

    add_executable(${TARGET} ${SOURCES} ${SOURCES_NS})
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
    message("Python3_LIBRARIES: ${Python3_LIBRARIES}")
    
    function(LOCATE_BOOST_PYTHON Boost_Python3_VERSION_TAG)
        foreach(TAG
                "${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}"
                "${Python3_VERSION_MAJOR}"
                "-py${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}"
                "-${Python3_VERSION_MAJOR}${Python3_VERSION_MINOR}"
                "${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}"
                "-${Python3_VERSION_MAJOR}.${Python3_VERSION_MINOR}")
            find_package(Boost QUIET COMPONENTS "python${TAG}")
            if(${Boost_FOUND})
                message("Boost_Python3_VERSION_TAG: ${TAG}")
                set(Boost_Python3_VERSION_TAG ${TAG} PARENT_SCOPE)
                break()
            endif()
        endforeach()
    endfunction(LOCATE_BOOST_PYTHON)

    find_package(Python3 COMPONENTS Development Interpreter)

    LOCATE_BOOST_PYTHON(Boost_Python3_VERSION_TAG)
    find_package(Boost COMPONENTS program_options filesystem "python${Boost_Python3_VERSION_TAG}" "numpy${Boost_Python3_VERSION_TAG}")
    message("Boost_LIBRARIES: ${Boost_LIBRARIES}")

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
