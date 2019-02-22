#ifndef __SPHSOLVER__
#define __SPHSOLVER__

#include <iostream>
#include <tuple>
#include <assert.h>
#include <CompactNSearch>
#include <chrono>
// #include <math>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "json.hpp"
#include "helpers.hpp"
#include "particleAttrib.hpp"

using namespace RealOps;
using json = nlohmann::json;
using Uint = unsigned int;
// <Uint,Uint> : First Uint is the set number, Second Uint is the particel number


class SPHSolver{
public:
	// const json& simData;
	json& simData;

	std::map<std::tuple<int,int,int>,bool> ghostMap;	
	std::map<std::string,ParticleAttributes*>& pData;
	std::vector<std::string> setNames;
	NeighborhoodSearch* nsearch;

	Real currentTime;
    SPHSolver(json& _simData, std::map<std::string,ParticleAttributes*>& pData, Real dx, Real dt, Real smoothingLength, Uint dims);
	void neighborSearch();
	void reloadAndNeighborSearch();
	void marchTime(Uint t);

	void initializeMass();


	void setDiffusiveTerm();
	void setKernels();
	void setEOS();
	void setThermalConductivityModel();
	void setHeatEquationDiscretization();
	void setTemperatureEnthalpyRelation();
	void setViscosityConstantFormulation();
	void setViscosityFormulation();
	void setBodyForce();
	void setPressureGradientFormulation();
	void setSensorParticles();
	void addFluidInletParticles(int t);

	void addGhostParticles(int t);
	void initGhostMap();
	void releaseGhostMap();
	void trimGhostParticles();


	// void calculateDensity();
	// void calculatePressureForce();
	// void calculateViscousForce();
	// void calculateHeatTransfer();
private:
	const Real3 zerovec{0.0,0.0,0.0};
	const Real3x3 zeromat{zerovec,zerovec,zerovec};
	const Real3x3x3 zeroijk{zeromat,zeromat,zeromat};
	const Real dt,dx,smoothingLength;
	const Uint dims;

	void setGhostParticleTemperatures();
	void setInitialConfigurationNeighbors();
	void setInitialDeformation();
	void computeInteractions(Uint t);
	void fixedPointIteration(Uint t);
	void XSPH(Uint t);
	void polarDecompose(MatrixXd& F, MatrixXd& R, MatrixXd& U);
	void getStretch(MatrixXd& U, Real3& lambdas, Real3x3& dirs, Real3x3& dirs_before);
	void smearDefGrad(Uint t);
	void smearStress(Uint t);
	void smearVelocity(Uint t);
	void computeDeformationGradient(Uint t);
	void computeDeformationGradient2(Uint t);

	std::function<Real(std::string type, Real T)> thermalConductivity;

	std::function<Real(Real3x3 L2_i, Real3 gradT_i,
					   Real Ti, Real Tj, Real mj, Real rho_i, Real rho_j, Real ki, Real kj,
					   Real3 relpos, Real3 reldir,
					   Real dist, Real3 gWij, Real vol_j)> heatTransfer;

	std::function<Real(Real)>															 TvsH;
	std::function<Real(Real,Real,Real,Real,Real,Real)>                					 EOS;
	std::function<Real(Real,Real)>         							  					 W_ij;
	std::function<Real3(Real,Real3,Real)>  							 					 gW_ij;
	std::function<Real3(Real,Real3,Real)>  							 			   	     gW_ij_nd;
	std::function<Real(Real,Real,Real,Real,Real,Real,Real,Real,Real3,Real3,Real3,Real3)> diffusiveTerm_ij;
	std::function<Real(Real,Real)>									  					 viscosityConstant;
	std::function<Real3()>                 							 			    	 bodyForceAcc_i;
	std::function<Real3(Real,Real3,Real3,Real,Real,Real,Real3,Real)>  					 viscosityAcc_ij;
	std::function<Real3(Real,Real3,Real,Real,Real,Real,Real)>  	  	  					 pressureAcc_ij;

	unsigned int totParticles;
	std::map<std::string,int> ids;



};

#endif
