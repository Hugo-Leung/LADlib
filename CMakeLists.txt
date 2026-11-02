# Software for LAD (Large Angle Detector) experiments in Hall C

cmake_minimum_required(VERSION 3.5)

project(LADlib LANGUAGES CXX)

if(NOT TARGET hcana:HallC)
  if(DEFINED ENV{HCANA})
    set(CMAKE_PREFIX_PATH "$ENV{HCANA}" "${CMAKE_PREFIX_PATH}")
  endif()

  find_package(hcana 0.95 REQUIRED)
endif()

include(PoddCMakeEnv)

set_diagnostic_flags(WALL WEXTRA)

add_subdirectory(src)
