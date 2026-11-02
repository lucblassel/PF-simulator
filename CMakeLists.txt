cmake_minimum_required(VERSION 3.20)
project(pf_alisim_native LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_NO_PACKAGE_REGISTRY ON)
# Cargo links these static libraries into a position-independent executable.
# Set this before creating any targets, including IQ-TREE's dependencies.
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Keep upstream's legacy build settings inside its own directory.
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)
set(BUILD_LIB ON CACHE BOOL "Build IQ-TREE as a library" FORCE)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Use static libraries" FORCE)
set(IQTREE_FLAGS "single" CACHE STRING "Sequential IQ-TREE build" FORCE)
set(USE_CMAPLE OFF CACHE BOOL "Disable unused inference dependency" FORCE)
set(USE_CMAPLE_AA OFF CACHE BOOL "Disable unused inference dependency" FORCE)
set(USE_LSD2 OFF CACHE BOOL "Disable dating" FORCE)
set(USE_TERRAPHAST OFF CACHE BOOL "Disable terraces" FORCE)
set(USE_BOOSTER OFF CACHE BOOL "Disable bootstrap support" FORCE)
add_subdirectory(dependencies/iqtree2 EXCLUDE_FROM_ALL)
include(native/memory_io.cmake)

# main.cpp supplies globals needed by the library as well as a CLI entry point.
# Rename that entry point so our executable can provide main(). No vendor edits.
target_compile_definitions(main PRIVATE main=pf_unused_iqtree_main)

# This revision has unguarded OpenMP calls in its sequential build.
foreach(iqtree_target simulator tree)
    target_compile_options(${iqtree_target} PRIVATE
        -include ${CMAKE_CURRENT_SOURCE_DIR}/native/iqtree_sequential_compat.hpp)
endforeach()

add_library(pf_alisim_bridge STATIC native/alisim_bridge.cpp)
get_directory_property(IQTREE_INCLUDES DIRECTORY dependencies/iqtree2 INCLUDE_DIRECTORIES)
get_directory_property(IQTREE_DEFINITIONS DIRECTORY dependencies/iqtree2 COMPILE_DEFINITIONS)
target_include_directories(pf_alisim_bridge PRIVATE ${IQTREE_INCLUDES})
target_compile_definitions(pf_alisim_bridge PRIVATE ${IQTREE_DEFINITIONS})
target_link_libraries(pf_alisim_bridge PRIVATE iqtree2)
# Cargo merges these transitive static archives into one native library.
function(collect_archives target)
    get_property(visited GLOBAL PROPERTY PF_ARCHIVES)
    if("${target}" IN_LIST visited)
        return()
    endif()
    get_target_property(kind ${target} TYPE)
    if(kind STREQUAL "STATIC_LIBRARY")
        set_property(GLOBAL APPEND PROPERTY PF_ARCHIVES ${target})
    endif()
    get_target_property(dependencies ${target} LINK_LIBRARIES)
    foreach(dependency IN LISTS dependencies)
        if(TARGET "${dependency}")
            collect_archives(${dependency})
        endif()
    endforeach()
endfunction()
collect_archives(pf_alisim_bridge)
get_property(archives GLOBAL PROPERTY PF_ARCHIVES)
set(manifest "")
foreach(archive IN LISTS archives)
    string(APPEND manifest "$<TARGET_FILE:${archive}>\n")
endforeach()
file(GENERATE OUTPUT "${CMAKE_BINARY_DIR}/native-archives.txt" CONTENT "${manifest}")
