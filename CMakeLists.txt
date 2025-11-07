cmake_minimum_required(VERSION 3.30)

project(Stat VERSION 1.0)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# --- Executable ---
add_executable(Stat)

# Sources
file(GLOB_RECURSE Stat_SOURCES CONFIGURE_DEPENDS 
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/*.inl"
)
target_sources(Stat PRIVATE 
	${Stat_SOURCES}
)

# Includes
target_include_directories(Stat PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/src
    include/
)

# --- Compiler warnings ---
if(MSVC)
    target_compile_options(Stat PRIVATE
        /W3
        /MP
        /WX
        /ZI
    )
else()
    target_compile_options(Stat PRIVATE 
        -Wall 
        -Wextra 
        -Wpedantic
    )
endif()

# --- Visual Studio niceties ---
if(MSVC)
	set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY VS_STARTUP_PROJECT Stat)
	source_group(TREE "${CMAKE_CURRENT_SOURCE_DIR}/src" PREFIX "Source Files" FILES ${Stat_SOURCES})
	set_property(GLOBAL PROPERTY USE_FOLDERS ON)
endif()