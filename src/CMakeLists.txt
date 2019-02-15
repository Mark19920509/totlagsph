
add_executable(FDMSPH
    main.cpp
    json.hpp
    sphsolver.hpp
    sphsolver.cpp
    sphsolver2.cpp
    particleAttrib.hpp
    particleAttrib.cpp
    helpers.hpp
    filewriter.hpp
    poissonsample.hpp
    geometrygeneration.cpp
    geometrygeneration.hpp
    vvector.h
)

target_link_libraries(FDMSPH CompactNSearch)
include_directories(${PROJECT_PATH}/extern/CompactNSearch/include)

find_package(hdf5 REQUIRED)
if(HDF5_INCLUDE_PATH)
    message("Found hdf5 at ${HDF5_INCLUDE_PATH}")
endif()
include_directories ( "${HDF5_INCLUDE_PATH}" )
target_link_libraries(FDMSPH "${HDF5_LIBRARY}")

find_package(H5Part REQUIRED)
include_directories ( "${H5PART_INCLUDE_PATH}" )
target_link_libraries(FDMSPH "${H5PART_LIBRARIES}" "-I/usr/local/hdf5/include -L/usr/local/hdf5/lib -lH5Part -lhdf5")
if(H5PART_LIBRARIES)
	message("Found H5Part at ${H5PART_INCLUDE_PATH}")
endif()


