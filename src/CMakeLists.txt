file(GLOB sources *.cpp)
file(GLOB headers *.h *.hpp)

list(REMOVE_ITEM sources main.cpp)

add_compile_options(-Wno-subobject-linkage) # Suppress bifrost warning

add_library(klue_core ${sources} ${headers})
target_include_directories(klue_core PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(klue main.cpp)

find_package( Threads REQUIRED )
ExternalProject_Get_Property(bifrost install_dir)

target_link_libraries(klue klue_core pthread ${install_dir}/build/src/libbifrost.a)

if(LINK MATCHES static)
    set(BUILD_SHARED_LIBS OFF)
    set(HDF5_USE_STATIC_LIBRARIES 1)

    if (UNIX AND NOT APPLE)
        #set(CMAKE_EXE_LINKER_FLAGS  "-static -static-libgcc -static-libstdc++")
        set(CMAKE_EXE_LINKER_FLAGS  "-static -static-libstdc++")
        SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
        set(CMAKE_EXE_LINKER_FLAGS  "-static -static-libgcc -static-libstdc++")
    endif(UNIX AND NOT APPLE)

    SET_TARGET_PROPERTIES(klue klue_core PROPERTIES LINK_SEARCH_END_STATIC 1)
endif(LINK MATCHES static)


if (NOT ZLIBNG)
find_package( ZLIB REQUIRED )

if ( ZLIB_FOUND )
    include_directories( ${ZLIB_INCLUDE_DIRS} )
    target_link_libraries(klue klue_core ${ZLIB_LIBRARIES})
else()
    message(FATAL_ERROR "zlib not found. Required for to output files" )
endif( ZLIB_FOUND )
endif(NOT ZLIBNG)

if (ZLIBNG)
    add_definitions(-DKLUE_USE_ZLIB_NG)
    include_directories(../ext/zlib-ng/zlib-ng)
    include_directories(../ext/zlib-ng)
    if(WIN32)
    target_link_libraries(klue ${CMAKE_CURRENT_SOURCE_DIR}/../ext/zlib-ng/zlib-ng/libz.lib)
    else()
    target_link_libraries(klue ${CMAKE_CURRENT_SOURCE_DIR}/../ext/zlib-ng/zlib-ng/libz.a)
    endif(WIN32)
endif(ZLIBNG)

if(LINK MATCHES static)
    if (UNIX AND NOT APPLE)
        target_link_libraries(klue librt.a)
    endif()
else()
    if (UNIX AND NOT APPLE)
        target_link_libraries(klue rt)
    endif()
endif(LINK MATCHES static)


install(TARGETS klue DESTINATION "${CMAKE_INSTALL_BINDIR}")
