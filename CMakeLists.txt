cmake_minimum_required(VERSION 3.22)
project(rtk LANGUAGES CXX)

add_library(rtk_core STATIC
    RTK/coordinate.cpp
    RTK/decode.cpp
    RTK/error_correction.cpp
    RTK/matrix.cpp
    RTK/obs.cpp
    RTK/rtk_ambiguity_resolution.cpp
    RTK/rtk_estimation.cpp
    RTK/rtk_mathematical_model.cpp
    RTK/rtk_prepare.cpp
    RTK/rtk_stochastic_model.cpp
    RTK/rtk_stream.cpp
    RTK/satpos.cpp
    RTK/spp.cpp
    RTK/station_satpos.cpp
    RTK/stream_decode.cpp
    RTK/timeTransform.cpp
)
target_include_directories(rtk_core PUBLIC RTK)
target_compile_features(rtk_core PUBLIC cxx_std_20)
set_target_properties(rtk_core PROPERTIES CXX_EXTENSIONS OFF)
if(WIN32)
    target_link_libraries(rtk_core PUBLIC ws2_32)
    target_compile_definitions(rtk_core PUBLIC WIN32_LEAN_AND_MEAN NOMINMAX)
endif()
if(MSVC)
    target_compile_options(rtk_core PUBLIC /utf-8)
    target_compile_definitions(rtk_core PUBLIC _CRT_SECURE_NO_WARNINGS)
endif()

# spp_main.cpp 是旧单站入口，不能与 rtk_main.cpp 同时链接。
add_executable(RTK RTK/rtk_main.cpp)
target_link_libraries(RTK PRIVATE rtk_core)

include(CTest)
if(BUILD_TESTING)
    add_executable(socket_probe portability_tests/socket_probe.cpp)
    target_link_libraries(socket_probe PRIVATE rtk_core)
    find_package(Python3 COMPONENTS Interpreter REQUIRED)
    add_test(NAME socket_transport
        COMMAND Python3::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/portability_tests/socket_transport.py
            $<TARGET_FILE:socket_probe>)
    if(UNIX)
        add_test(NAME stream_shutdown
            COMMAND Python3::Interpreter ${CMAKE_CURRENT_SOURCE_DIR}/portability_tests/stream_shutdown.py
                $<TARGET_FILE:RTK>)
    endif()
endif()
