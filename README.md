# fdmsph

# Consistent-SPH Implementation in CPU.


This code in based on the Laplacian Correction Matrix by Fatehi et al.
(http://www.sciencedirect.com/science/article/pii/S0898122110009004)

Sovles the issue in traditional SPH schemes, where each particle needs many neighbors to even provide 0th order consistency.

The popular linear algebra library `Eigen`, and the CompactNSearch algorithm by Dan Koschier was used (https://github.com/InteractiveComputerGraphics/CompactNSearch) extensively.

# Build instructions
Install libraries : 
``Eigen, HDF5 ,HDPART``
Clone the directory :
``git clone https://github.com/changyoonpark/fdmsph``
Create a build directory, use CMake to build.
``mkdir build && cd build``
``cmake .. && make``

For faster performance, you might want to use the `-Ofast` option with `ccmake .`.

# Drawback of the Method
- For each particle to be "complete" (the Laplacian Correction Matrix is invertible), there must exist at least one neighbor in each quadrant of the particle (at least 8 neighbors). Otherwise, the matirx is singular and the correction should not be used.


