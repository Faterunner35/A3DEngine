cmake_minimum_required(VERSION 3.20)
project(A3DEngine LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Release по умолчанию для одноконфигурационных генераторов (Ninja, Makefiles)
if(NOT CMAKE_CONFIGURATION_TYPES AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

option(A3D_NATIVE     "Оптимизировать под текущий процессор (-march=native, /arch:AVX2)" OFF)
option(A3D_SANITIZERS "AddressSanitizer + UBSan (GCC/Clang, Debug)" OFF)

# -----------------------------------------------------------------------------
#  Ядро: header-only математика, геометрия, ECS
# -----------------------------------------------------------------------------
add_library(a3d_core INTERFACE)
target_include_directories(a3d_core INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/SegmentTest)

if(MSVC)
    target_compile_options(a3d_core INTERFACE /utf-8 /W3 /permissive- /Zc:__cplusplus)
    target_compile_definitions(a3d_core INTERFACE _USE_MATH_DEFINES NOMINMAX)
    if(A3D_NATIVE)
        target_compile_options(a3d_core INTERFACE /arch:AVX2)
    endif()
else()
    target_compile_options(a3d_core INTERFACE -Wall -Wextra -Wno-unused-parameter)
    if(A3D_NATIVE)
        target_compile_options(a3d_core INTERFACE -march=native)
    endif()
    if(A3D_SANITIZERS)
        target_compile_options(a3d_core INTERFACE -fsanitize=address,undefined -fno-omit-frame-pointer)
        target_link_options(a3d_core INTERFACE -fsanitize=address,undefined)
    endif()
endif()

# -----------------------------------------------------------------------------
#  engine_tests — тесты, бенчмарк и стресс-нагрузка EntityPool + TransformComponent
# -----------------------------------------------------------------------------
add_executable(engine_tests SegmentTest/main.cpp)
target_link_libraries(engine_tests PRIVATE a3d_core)

# -----------------------------------------------------------------------------
#  CTest: только проверки корректности, без бенчмарка
# -----------------------------------------------------------------------------
enable_testing()
add_test(NAME transform_ecs COMMAND engine_tests --tests-only)

# Стартовый проект в Visual Studio
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT engine_tests)
