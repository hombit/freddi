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
    cpp/include/exceptions.hpp
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

if(NOT SKBUILD)
    CREATE_EXE("")
    CREATE_EXE("-ns")
endif()


function(CREATE_UNIT_TEST targ)
    set(TARGET test_${targ})

    find_package(Boost 1.57.0 COMPONENTS unit_test_framework program_options filesystem REQUIRED)

    add_executable(${TARGET} cpp/test/${targ}.cpp ${MIN_SRC} ${IO_SRC} ${NS_MIN_SRC} ${NS_IO_SRC})

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

if(NOT SKBUILD AND NOT STATIC_LINKING)
    CREATE_TESTS()
    enable_testing()
endif()


function(CREATE_PY_EXT)
    set(TARGET _freddi)

    find_package(PythonExtensions REQUIRED)
    find_package(NumPy REQUIRED)

    set(PYTHON_EXE "${PYTHON_PREFIX}/bin/python")
    set(PYTHON_VERSION_SCRIPT "import sys; vi = sys.version_info; print(f'{vi.major}{vi.minor}')")
    execute_process(
            COMMAND
            ${PYTHON_EXE} -c "${PYTHON_VERSION_SCRIPT}"
            OUTPUT_STRIP_TRAILING_WHITESPACE
            OUTPUT_VARIABLE PYTHON_VERSION
    )

    find_package(Boost REQUIRED COMPONENTS python${PYTHON_VERSION} numpy${PYTHON_VERSION} filesystem)

    add_library(${TARGET} MODULE ${MIN_SRC} ${NS_MIN_SRC} ${PYWRAP_SRC})
    target_include_directories(${TARGET} PUBLIC ${Boost_INCLUDE_DIRS} ${NumPy_INCLUDE_DIRS})
    target_link_libraries(${TARGET} ${Boost_LIBRARIES})
    python_extension_module(${TARGET})
    install(TARGETS ${TARGET} LIBRARY DESTINATION python/freddi)
    
endfunction(CREATE_PY_EXT)

if(SKBUILD)
    CREATE_PY_EXT()
endif()
