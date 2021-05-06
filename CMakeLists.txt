cmake_minimum_required(VERSION 3.16)
project(suloea)
set (CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(AEOLUS_DIR third_party/aeolus-0.9.9)

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