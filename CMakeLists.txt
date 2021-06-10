cmake_minimum_required(VERSION 3.1)
message("CMake version: ${CMAKE_VERSION}")
project(freddi)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

include_directories("${PROJECT_SOURCE_DIR}/cpp/include")

set(STATIC_LINKING FALSE CACHE BOOL "Build a static binary?")
if(STATIC_LINKING)
    set(Boost_USE_STATIC_LIBS "ON")
endif(STATIC_LINKING)

set(MIN_SRC
    cpp/src/arguments.cpp
    cpp/src/freddi_evolution.cpp
    cpp/src/freddi_state.cpp
    cpp/src/geometry.cpp
    cpp/src/nonlinear_diffusion.cpp
    cpp/src/opacity_related.cpp
    cpp/src/orbit.cpp
    cpp/src/passband.cpp
    cpp/src/rochelobe.cpp
    cpp/src/spectrum.cpp
    cpp/src/star.cpp
    cpp/src/util.cpp

    cpp/include/arguments.hpp
    cpp/include/constants.hpp
    cpp/include/freddi_evolution.hpp
    cpp/include/freddi_state.hpp
    cpp/include/geometry.hpp
    cpp/include/gsl_const_cgsm.h
    cpp/include/nonlinear_diffusion.hpp
    cpp/include/opacity_related.hpp
    cpp/include/orbit.hpp
    cpp/include/passband.hpp
    cpp/include/rochelobe.hpp
    cpp/include/spectrum.hpp
    cpp/include/star.hpp
    cpp/include/unit_transformation.hpp
    cpp/include/util.hpp
    )

set(IO_SRC
    cpp/src/options.cpp
    cpp/src/output.cpp
    cpp/include/options.hpp
    cpp/include/output.hpp
    )

set(NS_MIN_SRC
    cpp/src/ns/ns_arguments.cpp
    cpp/src/ns/ns_evolution.cpp
    cpp/include/ns/ns_arguments.hpp
    cpp/include/ns/ns_evolution.hpp
    )

set(NS_IO_SRC
    cpp/src/ns/ns_options.cpp
    cpp/src/ns/ns_output.cpp
    cpp/include/ns/ns_options.hpp
    cpp/include/ns/ns_output.hpp
    )

set(APP_SRC
    cpp/include/application.hpp
    )

set(PYWRAP_SRC
    cpp/pywrap/converters.cpp
    cpp/pywrap/converters.hpp
    cpp/pywrap/module.cpp
    cpp/pywrap/pywrap_arguments.cpp
    cpp/pywrap/pywrap_arguments.hpp
    cpp/pywrap/pywrap_freddi_evolution.cpp
    cpp/pywrap/pywrap_freddi_evolution.hpp
    cpp/pywrap/pywrap_freddi_state.cpp
    cpp/pywrap/pywrap_freddi_state.hpp
    cpp/pywrap/util.hpp
    )

set(SHARED_SRC ${MIN_SRC} ${IO_SRC} ${NS_MIN_SRC} ${NS_IO_SRC})


function(CREATE_SHARED)
    set(TARGET freddi-dll)
    find_package(Boost 1.57.0 COMPONENTS program_options filesystem REQUIRED)

    add_library(${TARGET} SHARED ${SHARED_SRC})
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

    find_package(Boost 1.57.0 COMPONENTS program_options filesystem REQUIRED)

    set(SOURCES ${MIN_SRC} ${IO_SRC} ${APP_SRC} cpp/main${targ}.cpp)
    if (targ STREQUAL "-ns")
        list(APPEND SOURCES ${NS_MIN_SRC} ${NS_IO_SRC})
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


function(CREATE_UNIT_TEST targ)
    set(TARGET test_${targ})

    find_package(Boost 1.57.0 COMPONENTS unit_test_framework program_options filesystem REQUIRED)

    add_executable(${TARGET} cpp/test/${targ}.cpp ${SHARED_SRC})

    target_include_directories(${TARGET} PUBLIC ${Boost_INCLUDE_DIR})
    target_link_libraries(${TARGET} ${Boost_LIBRARIES})

    add_test(${TARGET} ${TARGET})
endfunction()

function(CREATE_TESTS)
    file(GLOB TEST_SOURCES cpp/test/*.cpp)

    foreach(SOURCE_PATH IN LISTS TEST_SOURCES)
        get_filename_component(TEST_NAME ${SOURCE_PATH} NAME_WE)
        CREATE_UNIT_TEST(${TEST_NAME})
    endforeach()
endfunction()

CREATE_TESTS()
enable_testing()


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
    find_package(Boost COMPONENTS "python${Boost_Python3_VERSION_TAG}" "numpy${Boost_Python3_VERSION_TAG}")
    message("Boost_LIBRARIES: ${Boost_LIBRARIES}")

    if (${Python3_FOUND} AND ${Boost_FOUND})
        set(TARGET _freddi)
        Python3_add_library(${TARGET} MODULE ${MIN_SRC} ${NS_MIN_SRC} ${PYWRAP_SRC})
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
