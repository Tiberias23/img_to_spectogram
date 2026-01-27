cmake_minimum_required(VERSION 3.25)
project(img_to_spectrogram LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
find_package(OpenMP REQUIRED)

# ---------------------------------------------------------
# Hauptprogramm
# ---------------------------------------------------------
add_executable(img_to_spectrogram
        src/main.cpp
)

target_include_directories(img_to_spectrogram PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(img_to_spectrogram PRIVATE OpenMP::OpenMP_CXX)

if (CMAKE_BUILD_TYPE STREQUAL "Release")
    message("Setting optimization flags for Release build :3")
    target_compile_options(img_to_spectrogram PRIVATE $<$<CONFIG:Release>:-O3 -march=native -ffast-math -fopenmp>)
else ()
    message("a few optimisation flags set for Debug build for better comparison :3")
    # Optimierungen, aber Debug-Symbole behalten
    target_compile_options(img_to_spectrogram PRIVATE
            -O2                 # Optimierungen, aber nicht so aggressiv wie O3
            -g                  # Debug-Symbole
            -march=native
            -ffast-math
            -fopenmp
    )
endif()
