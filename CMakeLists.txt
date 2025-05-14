cmake_minimum_required(VERSION 3.10)
project(klue)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Silence Bifrost-specific warnings
add_compile_options(-Wno-subobject-linkage -Wno-deprecated-declarations)

# Include Bifrost headers (where Kmer.hpp lives)
include_directories(/usr/local/include/bifrost)

# If klue has a src/ folder, include that too
include_directories(${PROJECT_SOURCE_DIR}/src)

# Add your source files here
file(GLOB KLUE_SRC
    src/*.cpp
    src/*.h
)

add_library(klue_core STATIC ${KLUE_SRC})

# Link Bifrost
link_directories(/usr/local/lib)

add_executable(klue src/main.cpp)
target_link_libraries(klue klue_core bifrost)

find_package(ZLIB REQUIRED)

target_link_libraries(klue klue_core bifrost ZLIB::ZLIB)

if (BUILD_TESTING)
    add_subdirectory(${PROJECT_SOURCE_DIR}/ext/catch)
    include_directories(${CATCH_INCLUDE_DIR})
    add_subdirectory(unit_tests)
endif()

option(BUILD_FUNCTESTING "Build functional tests." OFF)

if (BUILD_FUNCTESTING)
    add_subdirectory(func_tests)
    message("Functional testing enabled.")
    add_custom_target(test
        COMMAND /bin/bash ./func_tests/runtests.sh
        DEPENDS ./src/klue
    )
endif()
