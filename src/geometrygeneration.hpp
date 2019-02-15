#include <CompactNSearch>

#include <iostream>
#include <fstream>
#include "json.hpp"
#include "helpers.hpp"
#include "filewriter.hpp"
#include "particleAttrib.hpp"
#include "poissonsample.hpp"
#include <cstdlib>
#include <random>  


namespace GeometryGeneration{

    void randomSphere(Real3 center, Real radius, Uint n, Real dx, std::vector<Real>& x, std::vector<Real>& y, std::vector<Real>& z);
    void randomBox(Real3 p1, Real3 p2, Uint n, Real dx, std::vector<Real>& x, std::vector<Real>& y, std::vector<Real>& z);
    void uniformCylinder(Real3 origin, Real length, Real radius, Uint nr, Uint n2, Real dx,
        std::vector<Real>& x,
        std::vector<Real>& y,
        std::vector<Real>& z);

    void randomCylinder(Real3 origin, Real length, Real radius, Uint n, Real dx, std::vector<Real>& x, std::vector<Real>& y, std::vector<Real>& z);

    void generateGeometry(json& simDataInput);
}