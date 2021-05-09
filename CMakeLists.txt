cmake_minimum_required(VERSION 3.16)
project(suloea)

# External configuration CMake scripts
set (CMAKE_MODULE_PATH "${CMAKE_MODULE_PATH};${CMAKE_CURRENT_SOURCE_DIR}/cmake")
include (SuloeaConfig)

set(AEOLUS_DIR third_party/aeolus-0.9.9)

### LV2
set (LV2PLUGIN_PRJ_NAME "${PROJECT_NAME}_lv2")

# Set the build directory as <build_dir>/lv2/<plugin_name>.lv2/
set (PROJECT_BINARY_DIR "${PROJECT_BINARY_DIR}/${PROJECT_NAME}.lv2")

# LV2 plugin specific settings
include (LV2Config)

add_library (${LV2PLUGIN_PRJ_NAME} MODULE 
    lv2/suloea.cpp
    ${AEOLUS_DIR}/source/addsynth.cc
    ${AEOLUS_DIR}/source/rankwave.cc
    ${AEOLUS_DIR}/source/division.cc
    ${AEOLUS_DIR}/source/asection.cc
    ${AEOLUS_DIR}/source/scales.cc
    ${AEOLUS_DIR}/source/reverb.cc
    ${AEOLUS_DIR}/source/rngen.cc
    ${AEOLUS_DIR}/source/exp2ap.cc
)
target_include_directories(${LV2PLUGIN_PRJ_NAME} PRIVATE ${AEOLUS_DIR}/source)
target_include_directories(${LV2PLUGIN_PRJ_NAME} PRIVATE third_party/)

# Explicitely strip all symbols on Linux but lv2_descriptor()
# MacOS linker does not support this apparently https://bugs.webkit.org/show_bug.cgi?id=144555
if (${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    file(COPY lv2/lv2.version DESTINATION ${CMAKE_BINARY_DIR})
    target_link_libraries(${LV2PLUGIN_PRJ_NAME} "-Wl,--version-script=lv2.version")
endif()
if (MINGW)
    set_target_properties (${LV2PLUGIN_PRJ_NAME} PROPERTIES LINK_FLAGS "-static")
endif()

# Remove the "lib" prefix, rename the target name and build it in the .lv build dir
# <build_dir>/lv2/<plugin_name>_lv2.<ext> to
# <build_dir>/lv2/<plugin_name>.lv2/<plugin_name>.<ext>
set_target_properties (${LV2PLUGIN_PRJ_NAME} PROPERTIES PREFIX "")
set_target_properties (${LV2PLUGIN_PRJ_NAME} PROPERTIES OUTPUT_NAME "${PROJECT_NAME}")
set_target_properties (${LV2PLUGIN_PRJ_NAME} PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

# Generate *.ttl files from *.in sources,
# create the destination directory if it doesn't exists and copy needed files
file (MAKE_DIRECTORY ${PROJECT_BINARY_DIR})
configure_file (lv2/manifest.ttl.in ${PROJECT_BINARY_DIR}/manifest.ttl)
configure_file (lv2/${PROJECT_NAME}.ttl.in ${PROJECT_BINARY_DIR}/${PROJECT_NAME}.ttl)
add_custom_command(TARGET ${LV2PLUGIN_PRJ_NAME} POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_directory
        "${PROJECT_SOURCE_DIR}/third_party/stops-0.3.0"
        "${PROJECT_BINARY_DIR}/stops"
)

### Aeolus

find_package(Freetype REQUIRED)
find_package(ALSA REQUIRED)
find_package(Threads REQUIRED)

add_library(aeolus_x11 SHARED 
    ${AEOLUS_DIR}/source/styles.cc
    ${AEOLUS_DIR}/source/mainwin.cc
    ${AEOLUS_DIR}/source/midiwin.cc
    ${AEOLUS_DIR}/source/audiowin.cc
    ${AEOLUS_DIR}/source/instrwin.cc
    ${AEOLUS_DIR}/source/editwin.cc
    ${AEOLUS_DIR}/source/midimatrix.cc
    ${AEOLUS_DIR}/source/multislider.cc
    ${AEOLUS_DIR}/source/functionwin.cc
    ${AEOLUS_DIR}/source/xiface.cc
    ${AEOLUS_DIR}/source/addsynth.cc
)
target_link_libraries(aeolus_x11 PRIVATE Xft X11 Freetype::Freetype clxclient)
target_compile_definitions(aeolus_x11 PRIVATE VERSION="\\"0.9.9\\"")
target_compile_definitions(aeolus_x11 PRIVATE LIBDIR=".")
set_target_properties(aeolus_x11 PROPERTIES PREFIX "")

add_executable(aeolus)
target_sources(aeolus PRIVATE 
    ${AEOLUS_DIR}/source/main.cc
    ${AEOLUS_DIR}/source/audio.cc
    ${AEOLUS_DIR}/source/model.cc
    ${AEOLUS_DIR}/source/slave.cc
    ${AEOLUS_DIR}/source/imidi.cc
    ${AEOLUS_DIR}/source/addsynth.cc
    ${AEOLUS_DIR}/source/scales.cc
    ${AEOLUS_DIR}/source/reverb.cc
    ${AEOLUS_DIR}/source/asection.cc
    ${AEOLUS_DIR}/source/division.cc
    ${AEOLUS_DIR}/source/rankwave.cc
    ${AEOLUS_DIR}/source/rngen.cc
    ${AEOLUS_DIR}/source/exp2ap.cc
    ${AEOLUS_DIR}/source/lfqueue.cc
)

target_link_libraries(aeolus PRIVATE zita-alsa-pcmi clthreads dl ALSA::ALSA Threads::Threads jack)

target_compile_definitions(aeolus PRIVATE VERSION="\\"0.9.9\\"")
target_compile_definitions(aeolus PRIVATE LIBDIR=".")
