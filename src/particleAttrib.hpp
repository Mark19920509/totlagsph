#ifndef __PARTICLEATTRIB__
#define __PARTICLEATTRIB__

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <climits>
#include <random>
#include <omp.h>
#include "helpers.hpp"
#include "csv.hpp"

using namespace CompactNSearch;
using namespace RealOps;

using Uint = unsigned int;

class ParticleAttributes{

public:
	unsigned int numParticles;
	unsigned int setID;
	Real dx, smoothingLength;
	const json& simDataIn, parDataIn;

	Real oneOverLambda, eta;

	Real defaultVolume, defaultMass;
	Real inletWidth;
	Real inletSpeed;
	Real3 inletVelocity;
	Real3 inletTangent;
	Real3 inletCenter;

	std::vector<std::vector<std::tuple<std::string,Uint>>> nMap;
	std::vector<Real>  vol;
	std::vector<Real>  originVol;
	std::vector<Real>  particleDensity;
	std::vector<Real3> particleDensityGrad;
	std::vector<Real> mass;
	std::vector<Real3> pos;
	std::vector<Real3> posBefore;
	std::vector<Real3> psi;
	std::vector<Real3> originPos;
	std::vector<Real3> perturb;
	std::vector<Real3> shift;

	std::vector<Real3> vel;
	std::vector<Real3> velBefore;

	std::vector<Real3> acc;
	std::vector<Real3> accBefore;
	std::vector<std::array<Real3,4>> vel_k;

	std::vector<Real3> force;
	std::vector<Real3> bodyForce;
	std::vector<Real3> densGrad;

	std::vector<Real3x3> dispGrad;
	std::vector<Real3x3> defoGrad;
	std::vector<Real3x3> defoGrad_thermal;
	std::vector<Real3x3> defoGrad_elastic;
	std::vector<Real3x3> defoGrad_withoutThermal;
	std::vector<Real3> stretch_total;
	std::vector<Real3> stretch_elastic;

	std::vector<Real3> stretch_plastic;
	std::vector<Real3> stretch_plastic_before;
	std::vector<Real3> stretch_plastic_dot;
	std::vector<Real3> stretch_plastic_dot_before;
	
	std::vector<Real3x3> stretch_dirs;
	std::vector<Real3x3> stretch_dirs_before;

	std::vector<Real3x3> secondPKStress;
	std::vector<Real3> normalVec;
	std::vector<Real>  curvature;
	std::vector<Real3> tempGrad;
	std::vector<Real3x3> L;
	std::vector<Real3x3> L_o;
	std::vector<Real> conditionNumber;
	std::vector<Real> conditionNumber2;
	std::vector<Real3x3> L2;

	std::vector<Real3x3> tau;
	std::vector<Real3x3> tauDot;
	std::vector<std::array<Real3x3,4>> tau_k;


	std::vector<Real3x3x3> tauGrad;
	std::vector<Real3x3> velGrad;

	std::vector<Real>  temp;
	std::vector<Real>  tempBefore;
	std::vector<Real>  tempDot;
	std::vector<Real>  tempDotBefore;

	std::vector<Real>  enthalpy;
	std::vector<Real>  enthalpydot;
	std::vector<std::array<Real,4>> enthalpy_k;

	std::vector<Real>  dens;
	std::vector<Real>  densdot;
	std::vector<std::array<Real,4>> dens_k;

	std::vector<bool>  isActive;
	std::vector<bool>  isSensor;
	std::vector<bool>  isGhost;
	std::vector<bool>  isFS;
	std::vector<bool>  isSolid;
	std::vector<bool>  isThermalDirichlet;
	std::vector<Real>  heatSensed;
	std::vector<Real3> forceSensed;
	std::vector<std::string> type;

	ParticleAttributes(const json& inputData, const json& particleData);

	void clearAll();
	void addParticlesToFluid();
	void addDefaultFluidParticleAtPosition(Real3& posToAdd);

//  This should only work for fluid particles
	inline Real getFirstLame(){        return parDataIn["FirstLame"];}
	inline Real getShearModulus(){       return parDataIn["shearModulus"];}
	inline Real getT0(){ 		  		 return parDataIn["T0"];}
	inline Real3 getv0(){                return Real3{parDataIn["v0"][0],parDataIn["v0"][1],parDataIn["v0"][2]};}
	inline Real getSpecificHeat(){ 		 return parDataIn["specificHeat"];}
	inline Real getRho0(){        		 return parDataIn["rho0"];}
	inline Real getMaterialDensity(){    return parDataIn["rho_material"];}
	inline Real getMu(){				 return parDataIn["mu"];}
	inline Real getSurfaceTensionCoeff(){ return parDataIn["surfaceTensionCoeff"];}
	inline Real getSoundSpeed() {		 return parDataIn["soundSpeed"];}
	inline Real getOneOverLambda(){      return 1./(Real)parDataIn["lambda"];}
	inline Real getEta(){                return (Real)parDataIn["eta"];}
	inline Real getThermalExpansion() {  return parDataIn["thermalExpansion"];}
	inline Real getViscosity(){          return parDataIn["viscosity"];}
	inline std::string getType(){        return parDataIn["type"];}
private:

	std::vector<Real> x,y,z,vx,vy,vz;

	const Real3     zeroVector{0,0,0};
	const Real3x3   zeromat{Real3{0,0,0},Real3{0,0,0},Real3{0,0,0}};
	const Real3x3x3 zeroijk{ Real3x3{Real3{0,0,0},Real3{0,0,0},Real3{0,0,0}},
							 Real3x3{Real3{0,0,0},Real3{0,0,0},Real3{0,0,0}},
							 Real3x3{Real3{0,0,0},Real3{0,0,0},Real3{0,0,0}}} ;

	void initParticles();
	void fluidInit();
	void boundaryInit();
	void ghostInit();
	void addDefaultFluidParticleAtPosition(Real3& posToAdd, Real temp);
	void readInitialPlacement(std::string fileName);
};

#define __PARTICLEATTRIBDEF_

#endif
