# Install script for directory: /Users/SJCY/Library/Mobile Documents/com~apple~CloudDocs/Research/FDMSPH/extern/CompactNSearch/src/ExternalProject_CompactNSearch

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/Users/SJCY/Library/Mobile Documents/com~apple~CloudDocs/Research/FDMSPH/extern/install/CompactNSearch")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include" TYPE FILE FILES
    "/Users/SJCY/Library/Mobile Documents/com~apple~CloudDocs/Research/FDMSPH/extern/CompactNSearch/src/ExternalProject_CompactNSearch/include/CompactNSearch"
    "/Users/SJCY/Library/Mobile Documents/com~apple~CloudDocs/Research/FDMSPH/extern/CompactNSearch/src/ExternalProject_CompactNSearch/include/Config.h"
    "/Users/SJCY/Library/Mobile Documents/com~apple~CloudDocs/Research/FDMSPH/extern/CompactNSearch/src/ExternalProject_CompactNSearch/include/CompactNSearch.h"
    "/Users/SJCY/Library/Mobile Documents/com~apple~CloudDocs/Research/FDMSPH/extern/CompactNSearch/src/ExternalProject_CompactNSearch/include/PointSet.h"
    "/Users/SJCY/Library/Mobile Documents/com~apple~CloudDocs/Research/FDMSPH/extern/CompactNSearch/src/ExternalProject_CompactNSearch/include/DataStructures.h"
    )
endif()

if("${CMAKE_INSTALL_COMPONENT}" STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/lib" TYPE STATIC_LIBRARY FILES "/Users/SJCY/Library/Mobile Documents/com~apple~CloudDocs/Research/FDMSPH/extern/CompactNSearch/src/ExternalProject_CompactNSearch-build/libCompactNSearch.a")
  if(EXISTS "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libCompactNSearch.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libCompactNSearch.a")
    execute_process(COMMAND "/Applications/Xcode.app/Contents/Developer/Toolchains/XcodeDefault.xctoolchain/usr/bin/ranlib" "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libCompactNSearch.a")
  endif()
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/Users/SJCY/Library/Mobile Documents/com~apple~CloudDocs/Research/FDMSPH/extern/CompactNSearch/src/ExternalProject_CompactNSearch-build/demo/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/Users/SJCY/Library/Mobile Documents/com~apple~CloudDocs/Research/FDMSPH/extern/CompactNSearch/src/ExternalProject_CompactNSearch-build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
