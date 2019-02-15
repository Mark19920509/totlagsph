FIND_LIBRARY(H5PART_LIBRARY NAMES H5Part )

FIND_PATH(H5PART_INCLUDE_PATH NAMES H5Part.h )

# Need to provide the *_LIBRARIES
SET(H5PART_LIBRARIES ${H5PART_LIBRARY})

# handle the QUIETLY and REQUIRED arguments and set H5PART_FOUND to TRUE if 
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(H5Part DEFAULT_MSG H5PART_LIBRARY H5PART_INCLUDE_PATH)