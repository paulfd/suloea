# Configuration for this plugin
# TODO: generate version from git
set (LV2PLUGIN_VERSION_MINOR   1)
set (LV2PLUGIN_VERSION_MICRO   0)
set (LV2PLUGIN_NAME            "suloea.lv2")
set (LV2PLUGIN_URI             "https://github.com/paulfd/suloea")
set (LV2PLUGIN_REPOSITORY      "https://github.com/paulfd/suloea")
set (LV2PLUGIN_AUTHOR          "Paul Ferrand")
set (LV2PLUGIN_EMAIL           "paul@ferrand.cc")
set (LV2PLUGIN_SPDX_LICENSE_ID "GPLv3")

if (MSVC)
    set (LV2PLUGIN_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/lv2" CACHE STRING
    "Install destination for LV2 bundle [default: ${CMAKE_INSTALL_PREFIX}/lv2}]")
else()
    set (LV2PLUGIN_INSTALL_DIR "${CMAKE_INSTALL_PREFIX}/lib/lv2" CACHE STRING
    "Install destination for LV2 bundle [default: ${CMAKE_INSTALL_PREFIX}/lib/lv2}]")
endif()
