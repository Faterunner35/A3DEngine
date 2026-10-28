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

# Папки по группам в Solution Explorer (CMakePredefinedTargets и т. п.)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

# -----------------------------------------------------------------------------
#  Исходники движка: ВСЁ, что лежит в SegmentTest/ (включая подпапки).
#
#  Список файлов здесь не ведётся. Новый .cpp/.h/.hpp, добавленный в SegmentTest/
#  через Visual Studio, Проводник или git, попадает в сборку сам: CONFIGURE_DEPENDS
#  заставляет CMake перепроверять папку при каждой сборке и перегенерировать
#  проект, если состав файлов изменился.
#
#  Файлы с функцией main() — точки входа программ — перечислены в A3D_ENTRY_POINTS
#  и в библиотеку не идут. Новую программу добавлять туда же и ниже — add_executable.
# -----------------------------------------------------------------------------
set(A3D_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SegmentTest)

file(GLOB_RECURSE A3D_SOURCES CONFIGURE_DEPENDS
    ${A3D_SRC_DIR}/*.cpp
    ${A3D_SRC_DIR}/*.h
    ${A3D_SRC_DIR}/*.hpp
    ${A3D_SRC_DIR}/*.inl)

set(A3D_ENTRY_POINTS
    ${A3D_SRC_DIR}/main.cpp       # engine_tests
    ${A3D_SRC_DIR}/PgaDemo.cpp)   # pga_demo
list(REMOVE_ITEM A3D_SOURCES ${A3D_ENTRY_POINTS})

# -----------------------------------------------------------------------------
#  a3d_core — статическая библиотека движка. Заголовки входят в цель, чтобы
#  Visual Studio показывала их в дереве проекта и давала по ним IntelliSense.
# -----------------------------------------------------------------------------
add_library(a3d_core STATIC ${A3D_SOURCES})
source_group(TREE ${A3D_SRC_DIR} PREFIX "Engine" FILES ${A3D_SOURCES})

# Папку SegmentTest НЕ добавляем в пути поиска: в ней лежит свой Math.h, и на
# Windows (регистр не различается) <cmath> подхватил бы его вместо системного
# <math.h>. Заголовки подключаются в кавычках относительно своего файла.

if(MSVC)
    target_compile_options(a3d_core PUBLIC /utf-8 /W3 /permissive- /Zc:__cplusplus /MP)
    target_compile_definitions(a3d_core PUBLIC _USE_MATH_DEFINES NOMINMAX)
    if(A3D_NATIVE)
        target_compile_options(a3d_core PUBLIC /arch:AVX2)
    endif()
else()
    target_compile_options(a3d_core PUBLIC -Wall -Wextra -Wno-unused-parameter)
    if(A3D_NATIVE)
        target_compile_options(a3d_core PUBLIC -march=native)
    endif()
    if(A3D_SANITIZERS)
        target_compile_options(a3d_core PUBLIC -fsanitize=address,undefined -fno-omit-frame-pointer)
        target_link_options(a3d_core PUBLIC -fsanitize=address,undefined)
    endif()
endif()

# -----------------------------------------------------------------------------
#  Программы
# -----------------------------------------------------------------------------
# engine_tests — тесты, бенчмарк и стресс-нагрузка EntityPool + TransformComponent
add_executable(engine_tests ${A3D_SRC_DIR}/main.cpp)
target_link_libraries(engine_tests PRIVATE a3d_core)

# pga_demo — демонстрация PGA Cl(3,0,1) (бывший SegmentTest.cpp)
add_executable(pga_demo ${A3D_SRC_DIR}/PgaDemo.cpp)
target_link_libraries(pga_demo PRIVATE a3d_core)

source_group("Apps" FILES ${A3D_ENTRY_POINTS})

# -----------------------------------------------------------------------------
#  CTest: только проверки корректности, без бенчмарка
# -----------------------------------------------------------------------------
enable_testing()
add_test(NAME transform_ecs COMMAND engine_tests --tests-only)
add_test(NAME pga_demo      COMMAND pga_demo)

# Стартовый проект в Visual Studio
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT engine_tests)
