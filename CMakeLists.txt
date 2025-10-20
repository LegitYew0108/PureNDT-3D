project(PureNDT3D)
cmake_minimum_required(VERSION 3.10)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_library(PureNDT3D SHARED
	src/ndt_core.cpp
	src/ndt_calculator.cpp
	src/voxel.cpp
	src/ndt_optimizer.cpp
	src/ndt_matcher.cpp
)

target_include_directories(PureNDT3D PUBLIC
	include)

if(BUILD_TESTING)
    # Enable testing, and if so, download and configure GTest
    enable_testing()

    include(FetchContent)
    FetchContent_Declare(
      googletest
      URL https://github.com/google/googletest/archive/refs/tags/v1.14.0.zip
    )
    # For Windows: Prevent overriding the parent project's compiler/linker settings
    set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
    FetchContent_MakeAvailable(googletest)

    add_executable(interface_test
        test/interface_test.cpp)
    target_link_libraries(interface_test
        PureNDT3D
        GTest::gtest_main)

    include(GoogleTest)
    gtest_discover_tests(interface_test)
endif()
