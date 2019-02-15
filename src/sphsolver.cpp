#include "sphsolver.hpp"
#define CALC_HEAT 1

SPHSolver::SPHSolver(json& _simData, std::map<std::string,ParticleAttributes*>& _pData, Real _dx, Real _dt, Real _smoothingLength, Uint _dims)
// SPHSolver::SPHSolver(const json& _simData, std::map<std::string,ParticleAttributes*>& _pData)
: simData(_simData), pData(_pData), dx(_dx), dt(_dt), smoothingLength(_smoothingLength), dims(_dims)
{

	// Generate Particle Atrribution Arrays for the fluid
	nsearch = new NeighborhoodSearch((double)simData["smoothingLength"]*1.1,true);
	// load the particlesets onto the nsearcher

	std::cout << " ***** Initializing SPH Solver *****" << std::endl;
	currentTime = 0.0;

	for (const auto& pDatEntry : pData){
		ids[pDatEntry.first] = nsearch->add_point_set( pDatEntry.second->pos.front().data(), pDatEntry.second->pos.size(), true, true);
		std::cout << "... Particle set \"" << pDatEntry.first << "\" with " << pDatEntry.second->numParticles << " Particles Set Loaded onto CompactNSearch." << std::endl;
		if (pDatEntry.second->pos.size() != pDatEntry.second->numParticles)	assert(!(pDatEntry.second->pos.size() == pDatEntry.second->numParticles));
		setNames.push_back(pDatEntry.first);
	}


	totParticles = 0;
	Uint numSets = 0;

	for (const auto& pSet : nsearch->point_sets()){
		totParticles += pSet.n_points();
		numSets += 1;
	}

	std::cout << "... Loaded " << totParticles << " particles, with " << numSets << " point sets." << std::endl;

	// Set the SPH model.
	std::cout << "... defining models for SPH" << std::endl;
	setEOS();
	setDiffusiveTerm();
	setKernels();
	setBodyForce();
	setPressureGradientFormulation();
	setViscosityFormulation();
	setViscosityConstantFormulation();
	setEOS();
	setThermalConductivityModel();
	setHeatEquationDiscretization();
	setTemperatureEnthalpyRelation();

	//Set the sensor particles.
	setSensorParticles();

	// Initialize the mass / volume
	// neighborSearch();
	reloadAndNeighborSearch();
	setInitialConfigurationNeighbors();
	if( simData["useVariableVolume"] == "Yes" ) initializeMass();


	std::cout << "--- total number of particle sets : " << nsearch->n_point_sets() << std::endl;
	std::cout << "--- total number of particles     : " << totParticles << std::endl;

	// setInitialDeformation();
}

void SPHSolver::addFluidInletParticles(int t){

	if ( t % (int)( ((Real) (pData["fluid"]->dx)) / ((Real) (pData["fluid"]->inletSpeed) * (Real) simData["dt"] ) ) != 0 ) return;

	std::cout << "... Adding inlet particles" << std::endl;
	pData["fluid"]->addParticlesToFluid();

	std::cout << "... Freeing nsearch instance" << std::endl;
	delete nsearch;
	nsearch = new NeighborhoodSearch((double)simData["smoothingLength"],true);

	for (const auto& pDatEntry : pData){
		ids[pDatEntry.first] = nsearch->add_point_set( pDatEntry.second->pos.front().data(), pDatEntry.second->pos.size(), true, true);
		std::cout << "... Particle set \"" << pDatEntry.first << "\" with " << pDatEntry.second->numParticles << " Particles Set Loaded onto CompactNSearch." << std::endl;
		if (pDatEntry.second->pos.size() != pDatEntry.second->numParticles)	assert(!(pDatEntry.second->pos.size() == pDatEntry.second->numParticles));
	}

	// There is a bug in the nsearch code that freezes the code at this part.
	// nsearch->resize_point_set(ids["fluid"], 
	// 							  pData["fluid"]->pos.front().data(),
	// 							  pData["fluid"]->pos.size());
								
	// if (pData["fluid"]->pos.size() != pData["fluid"]->numParticles)	assert(!(pData["fluid"]->pos.size() == pData["fluid"]->numParticles));

	std::cout << "... Resized searching algorithm for the inlet particles" << std::endl;

}


void SPHSolver::initializeMass(){

	std::cout << "--- Initializing particle volume / mass." << std::endl;
	const Real smoothingLength = (Real) simData["smoothingLength"];
	Real totVol = 0.0;

	for (const auto& setName_i : setNames){

		const int setID_i = ids[setName_i];
		const auto& ps_i = nsearch->point_set(setID_i);

		// Real totmass = 0;
		for (int i = 0; i < ps_i.n_points(); ++i){
			Real kernelSum = 0.;
			for (const auto& setName_j : setNames){

				const int setID_j = ids[setName_j];
				const auto& ps_j = nsearch->point_set(setID_j);

				for (int _j = 0; _j < ps_i.n_neighbors(setID_j,i); _j++){
					Uint const j = ps_i.neighbor(setID_j, i, _j);

					Real3  relpos = sub(pData[setName_i]->pos[i],pData[setName_j]->pos[j]);
					Real     dist = length(relpos);

					if ( dist > smoothingLength ) continue;

					Real3  relvel = sub(pData[setName_i]->vel[i],pData[setName_j]->vel[j]);
					Real3  reldir = divide(relpos,dist);

					Real      Wij =  W_ij(dist, smoothingLength);
					kernelSum += Wij;
				}

			}
			kernelSum += W_ij(0.0,smoothingLength);

			pData[setName_i]->vol[i]  = (1.0/kernelSum);
			totVol += pData[setName_i]->vol[i];
			pData[setName_i]->mass[i] =  pData[setName_i]->dens[i] * pData[setName_i]->vol[i];

		}
  }
  simData["totalVolume"] = (Real)totVol;  

}

void SPHSolver::setInitialConfigurationNeighbors(){
	std::cout << "--- Initializing Initial Configuration Neighbor Map" << std::endl;
	const Real smoothingLength = (Real) simData["smoothingLength"];
	Real totVol = 0.0;

	for (const auto& setName_i : setNames){
		// std::cout << setName_i << std::endl;

		const int setID_i = ids[setName_i];
		if (setName_i != "fluid") continue;
		const auto& ps_i = nsearch->point_set(setID_i);

		for (int i = 0; i < ps_i.n_points(); ++i){
			pData[setName_i]->particleDensity[i] = 0;

			for (const auto& setName_j : setNames){

				if (setName_j != "fluid") continue;
				const int setID_j = ids[setName_j];
				const auto& ps_j = nsearch->point_set(setID_j);

				for (int _j = 0; _j < ps_i.n_neighbors(setID_j,i); _j++){
					Uint const j = ps_i.neighbor(setID_j, i, _j);					
					Real3  relpos = sub(pData[setName_i]->pos[i],pData[setName_j]->pos[j]);
					Real     dist = length(relpos);	
					Real      Wij = W_ij(dist, smoothingLength);					
					pData[setName_i]->nMap[i].push_back(std::tuple<std::string,Uint>{setName_j,j});
					if (dist < smoothingLength){
						pData[setName_i]->particleDensity[i] += pData[setName_j]->vol[j] * Wij; 
					}
				}

			}	

			pData[setName_i]->particleDensity[i] += pData[setName_i]->vol[i] * W_ij(0, smoothingLength);; 
		}
	}
}

void SPHSolver::setInitialDeformation(){
	for (const auto& setName_i : setNames){

		const int setID_i = ids[setName_i];
		const auto& ps_i = nsearch->point_set(setID_i);

		for (int i = 0; i < ps_i.n_points(); ++i){
				if(setName_i == "fluid" && pData[setName_i]->isSolid[i]){
					pData[setName_i]->pos[i] = add(pData[setName_i]-> originPos[i] , Real3{((0) * pData[setName_i]->originPos[i][0]  + 0.0  * pData[setName_i]->originPos[i][1]),
																							(0.01 * pData[setName_i]->originPos[i][0]   + 0.0  * pData[setName_i]->originPos[i][1]),0} );				
				}
		}

	}
}


void SPHSolver::setSensorParticles(){
	for(auto& sensor : simData["sensors"]){
		if ( sensor["type"] == "sensorZPlane" ){
			Real sensorLocation = (Real)sensor["coord"];
			std::cout << "--- Creating Sensor Plane (Z = " << sensor["coord"] << ")" << std::endl;
			for (const auto& setName : setNames){
				const int setID = ids[setName];
				const auto& ps = nsearch->point_set(setID);
				if (setName == "boundary"){
					#pragma omp parallel for num_threads(NUMTHREADS)
					for (int i = 0; i < ps.n_points(); ++i){
						if ( abs(pData[setName]->pos[i][2] - sensorLocation) < simData["smoothingLength"] ){
							pData[setName]->isSensor[i] = true;
						}
					}
				}
			}

	}}
}

void SPHSolver::setThermalConductivityModel(){
	if( simData["thermalConductivity"] == "Constant"){
		std::cout << "--- Thermal Conductivity Model : Constant conductivity" << std::endl;
		thermalConductivity = constantConductivity;
	} else{
		assert(1 &&& "!!! Unimplemented Thermal Conductivity Model.");
	}
}

void SPHSolver::setHeatEquationDiscretization(){

	if (simData["heatEqDiscretization"] == "Cleary"){
		std::cout << "--- Heat Equation Discretization Model : Monaghan-Cleary Formulation" << std::endl;
		// heatTransfer = cleary;
		assert(0 && "!!! Unimplemented Heat Equation Discretization.");
	} else if (simData["heatEqDiscretization"] == "Consistent"){
		std::cout << "--- Heat Equation Discretization Model : Consistent Renormalized Laplacian Formulation" << std::endl;
		heatTransfer = consistentHeatTransfer;
	} else{
		assert(1 && "!!! Unimplemented Heat Equation Discretization.");
	}
}


void SPHSolver::setTemperatureEnthalpyRelation(){
	if (simData["temperatureEnthalpyRelation"] == "T=H"){
		std::cout << "--- Temperature-Enthalpy Relation : T=H (Debugging ONLY)" << std::endl;
		TvsH = equal;
	} else{
		assert(1 && "!!! Unimplemented Temperature / Enthalpy Relation.");
	}
}

void SPHSolver::setEOS(){
	if (simData["EOS"] == "Linear"){
		std::cout << "--- EOS : Linear (P-V-T Variant)" << std::endl;
		EOS = linearEOS;
	} else if (simData["EOS"] == "Tait"){
		std::cout << "--- EOS : Tait (P-V Variant)" << std::endl;
		EOS = taitEOS;
	} else{
		assert(1 && "!!! Unimplemented Equation of state.");
	}
}


void SPHSolver::setDiffusiveTerm(){
	if (simData["diffusionModel"] == "DeltaSPH"){
		std::cout << "--- Diffusion Model : DeltaSPH" << std::endl;
		diffusiveTerm_ij = delta_SPH;
	} else if (simData["diffusionModel"] == "None"){
		std::cout << "--- Diffusion Model : None" << std::endl;
		diffusiveTerm_ij = [](Real  rho_i,      Real rho_j,
 						  	  Real  vol_j,      Real delta,
 						  	  Real  soundSpeed, Real smoothingLength,
 						  	  Real  mass_j,     Real  dist,
 						  	  Real3 relpos,    	Real3 gWij,
 						  	  Real3 densGrad_i, Real3 densGrad_j){return 0;};
	} else{
		assert(0 && "!!! Unimplemented Diffusion Model.");
	}
}

void SPHSolver::setKernels(){
	if (simData["kernelType"] == "Wendland"){
		if (simData["dimensions"] == 2){
			std::cout << "--- Kernel Type : Wendland (2D)" << std::endl;
			W_ij  =  W_Wendland_2D;
			gW_ij = gW_Wendland_2D;
		} else if (simData["dimensions"] == 3){
			std::cout << "--- Kernel Type : Wendland (3D)" << std::endl;
			W_ij  =  W_Wendland_3D_2h;
			gW_ij = gW_Wendland_3D_2h;
		}
	} 
	else if (simData["kernelType"] == "Wendland6"){
		if (simData["dimensions"] == 2){
			std::cout << "--- Kernel Type : Wendland 6 (2D)" << std::endl;
			W_ij  =  W_Wendland6_2D_h;
			gW_ij = gW_Wendland6_2D_h;
		} else{
			assert(0 && "!!! Unimplemented Kernel Type.");
		}
	}
	else if (simData["kernelType"] == "Wendland_2h"){
		if (simData["dimensions"] == 2){
			std::cout << "--- Kernel Type : Wendland (2D, 2h version)" << std::endl;
			W_ij  =  W_Wendland_2D_2h;
			gW_ij = gW_Wendland_2D_2h;
			gW_ij_nd = gW_Wendland_2D_2h_Nondim; 			
		} else if (simData["dimensions"] == 3){
			std::cout << "--- Kernel Type : Wendland (3D, 2h version)" << std::endl;
			W_ij  =  W_Wendland_3D_2h;
			gW_ij = gW_Wendland_3D_2h;
			gW_ij_nd = gW_Wendland_3D_2h_Nondim; 			
		} else{
			std::cout << "--- Kernel Type : Wendland (1D, 2h version)" << std::endl;
			W_ij  =  W_Wendland_1D_2h;
			gW_ij = gW_Wendland_1D_2h;
			gW_ij_nd = gW_Wendland_1D_2h_Nondim; 						
		}
	} else if (simData["kernelType"] == "Quintic_Spline"){
		if (simData["dimensions"] == 3){
			std:: cout << "--- Kernel Type : Quintic Spline (3D, h version)" << std::endl;
			W_ij  = W_QuinticSpline_3D;
			gW_ij = gW_QuinticSpline_3D;
		}
	} 
	else{
		assert(0 && "!!! Unimplemented Kernel Type.");
	}
}

void SPHSolver::setBodyForce(){
	if (simData["bodyForce"] == "Gravity"){
		std::cout << "--- Body Force Type : Gravity" << std::endl;
		bodyForceAcc_i = gravity_acc;
	} else if (simData["bodyForce"] == "None"){
		std::cout << "--- Body Force Type : None" << std::endl;
		bodyForceAcc_i = [](){return Real3{0.0,0.0,0.0};};
	} else{
		assert(1 && "!!! Unimplemented BodyForce Type.");
	}
}

void SPHSolver::setViscosityConstantFormulation(){
	if (simData["viscosityConstantFormulation"] == "temperatureDependent"){
	} else if (simData["viscosityConstantFormulation"] == "Fixed"){
		std::cout << "--- Viscosity Constant Formulation : Fixed Viscosity" << std::endl;
		viscosityConstant = fixedViscosity;
	}
}

void SPHSolver::setViscosityFormulation(){
	if (simData["viscosityFormulation"] == "Shao"){
		std::cout << "--- Viscosity Formulation : Shao" << std::endl;
		viscosityAcc_ij = viscosity_acc_ij_Shao;
	} else{
		assert(1 && "!!! Unimplemented viscosity formualtion.");
	}
}

void SPHSolver::setPressureGradientFormulation(){

	if ((int)simData["pressureFormulation"] == 1){
		std::cout << "--- Pressure Gradient Formulation : Type 1" << std::endl;
		pressureAcc_ij = pressure_acc_ij_1;
	} else if ((int)simData["pressureFormulation"] == 2){
		std::cout << "--- Pressure Gradient Formulation : Type 2" << std::endl;
		pressureAcc_ij = pressure_acc_ij_2;
	} else if ((int)simData["pressureFormulation"] == 3){
		std::cout << "--- Pressure Gradient Formulation : Type 3" << std::endl;
		pressureAcc_ij = pressure_acc_ij_3;
	} else{
		assert(1 && "!!! Unimplemented pressure gradient formualtion.");
	}

}

void SPHSolver::neighborSearch(){

	std::cout << "... Searching Neighbors" << std::endl;



	auto t0 = std::chrono::high_resolution_clock::now();

	// for (const auto& setName_i : setNames){
		
	// 			const int setID_i = ids[setName_i];
	// 			const auto& ps_i = nsearch->point_set(setID_i);
		
	// 			for (int i = 0; i < ps_i.n_points(); ++i){		
	// 				Real3 pos_i  = pData[setName_i]->pos[i];
	// 				std::cout << pos_i[0] << ", " << pos_i[1] << ", " << pos_i[2] << std::endl;
	// 			}
	// }

	nsearch->find_neighbors();
	std::cout << "--- Neighborhood search took " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0).count() << "ms" << std::endl;


	// I have no idea what this does. Documentation is unclear.
	// nsearch->z_sort();
	// for (const auto& pDatEntry : pData){
	// 	auto& ps = nsearch->point_set(ids[pDatEntry.first]);
	// 	ps.sort_field(pDatEntry.second-> pos.data());
	// 	ps.sort_field(pDatEntry.second-> vel.data());
	// 	ps.sort_field(pDatEntry.second-> acc.data());
	// 	ps.sort_field(pDatEntry.second->dens.data());
	// 	ps.sort_field(pDatEntry.second->temp.data());
	// 	ps.sort_field(pDatEntry.second->type.data());
	// }

}



void SPHSolver::polarDecompose(	MatrixXd& F, MatrixXd& R, MatrixXd& U){
	MatrixXd W = MatrixXd::Zero(3,3), SIGM = MatrixXd::Zero(3,3), V = MatrixXd::Zero(3,3);
	JacobiSVD<MatrixXd> svd_F(F, ComputeFullU | ComputeFullV);
	SIGM = svd_F.singularValues().asDiagonal();
	W    = svd_F.matrixU();
	V    = svd_F.matrixV();
	R    = W * V.transpose();
	U	 = V * SIGM * V.transpose();
}

void SPHSolver::getStretch(MatrixXd& U, Real3& lambdas, Real3x3& dirs, Real3x3& dirs_before){

	EigenSolver<MatrixXd> es(U);

	Real3x3 basis{zeromat};
	Real3   _lambdas{zerovec};

	toReal3(es.eigenvalues().real(), _lambdas);
	toReal3(es.eigenvectors().col(0).real(), basis[0]);
	toReal3(es.eigenvectors().col(1).real(), basis[1]);
	toReal3(es.eigenvectors().col(2).real(), basis[2]);

	//ugly stupid routine to sort the eigvecs
	Real a = std::abs(dot(basis[0], dirs_before[0])), b = std::abs(dot(basis[0], dirs_before[1])), c = std::abs(dot(basis[0], dirs_before[2]));
	if( a > b && a > c){
		// now[0] = before[0]
		dirs[0] = basis[0];
		lambdas[0] = _lambdas[0];
		Real d = std::abs(dot(basis[1], dirs_before[1])), e = std::abs(dot(basis[1], dirs_before[2]));
		if ( d > e ){
			// now[1] = before[1], now[2] = before[2]
			dirs[1] = basis[1]; dirs[2] = basis[2];
			lambdas[1] = _lambdas[1]; lambdas[2] = _lambdas[2];
		} else{
			// now[1] = before[2], now[2] = before[1]
			dirs[2] = basis[1]; dirs[1] = basis[2];
			lambdas[2] = _lambdas[1]; lambdas[1] = _lambdas[2];
		}
	} else if ( b > a && b > c){
		// now[0] = before[1]
		dirs[1] = basis[0];
		lambdas[1] = _lambdas[0];
		Real d = std::abs(dot(basis[1], dirs_before[0])), e = std::abs(dot(basis[1], dirs_before[2]));
		if ( d > e ){
			// now[1] = before[0], now[2] = before[2]
			dirs[0] = basis[1]; dirs[2] = basis[2];
			lambdas[0] = _lambdas[1]; lambdas[2] = _lambdas[2];
		} else{
			// now[1] = before[2], now[2] = before[0]
			dirs[2] = basis[1]; dirs[0] = basis[2];
			lambdas[2] = _lambdas[1]; lambdas[0] = _lambdas[2];
		}		
	} else{
		// now[0] = before[2]
		dirs[2] = basis[0];
		lambdas[2] = _lambdas[0];
		Real d = std::abs(dot(basis[1], dirs_before[0])), e = std::abs(dot(basis[1], dirs_before[1]));
		if ( d > e ){
			// now[1] = before[0], now[2] = before[1]
			dirs[0] = basis[1]; dirs[1] = basis[2];
			lambdas[0] = _lambdas[1]; lambdas[1] = _lambdas[2];
		} else{
			// now[1] = before[1], now[2] = before[0]
			dirs[1] = basis[1]; dirs[0] = basis[2];
			lambdas[1] = _lambdas[1]; lambdas[0] = _lambdas[2];
		}		
	}
}


void SPHSolver::computeDeformationGradient(Uint t){
	using namespace RealOps;
	const Real dx = (Real) simData["dx"];
	const Real dt = (Real) simData["dt"];
	const Uint dims = simData["dimensions"];
	const Real smoothingLength = (Real) simData["smoothingLength"];
	
	for (const auto& setName_i : setNames){


		const int setID_i = ids[setName_i];
		const auto& ps_i = nsearch->point_set(setID_i);
		
		// #pragma omp parallel for num_threads(NUMTHREADS)
		for (int i = 0; i < ps_i.n_points(); ++i){


			// This should only be computed for a solidified fluid particle.
			if ( !pData[setName_i]->isSolid[i] ) continue;		


			Real3& pos_i = pData[setName_i]->pos[i];
			// Renormalization Tensor (First Derivative)
			// This can claimed to be "inconsistent", but it is something that we accept
			MatrixXd _L_i(3,3); _L_i = MatrixXd::Zero(3,3);
			Real3x3& L_i = pData[setName_i]->L[i];
			L_i = zeromat;

			// Initialize Deformation Gradient and the cauchy stress tensor
			// if( t <= 100){
				// pData[setName_i]->temp[i] = t;
			// } 
			// else{
			// 	pData[setName_i]->temp[i] = std::max(
			// 								(Real)(10.0 * pData[setName_i]->originPos[i][1] * pData[setName_i]->originPos[i][1] * 1000.0 * 0.2 - 
			// 								10.0 * pData[setName_i]->originPos[i][1] * pData[setName_i]->originPos[i][1] * ((Real)t-1000.0) * 0.2),
			// 								(Real)0.0);
			// }

			Real thermalDeform = 1.0 + pData[setName_i]->getThermalExpansion() * pData[setName_i]->temp[i];
			pData[setName_i]->defoGrad[i] = zeromat;
			if ( pData[setName_i]->isSolid[i] && setName_i == "fluid" ){
				// pData[setName_i]->defoGrad_thermal[i] = mult((thermalDeform),identity());
				pData[setName_i]->defoGrad_thermal[i] = mult( thermalDeform ,identity());
			}

			// Compute the Renormalization Matrix / density gradient using the renormalization matrix.
			for (const auto& setName_j : setNames){

				const int setID_j = ids[setName_j];
				const auto& ps_j = nsearch->point_set(setID_j);

				for (int _j = 0; _j < ps_i.n_neighbors(setID_j,i); _j++){

					// Index of neighbor particle in setID_j
					Uint const j = ps_i.neighbor(setID_j, i, _j);

					// This should only be computed for a solidified fluid particle.
					if( !pData[setName_j]->isSolid[j] ) continue;

					// Define all the ingedients for the particle interaction
					Real3  relpos = sub(pData[setName_i]->pos[i],pData[setName_j]->pos[j]);
					Real     dist = length(relpos);

					if ( dist > smoothingLength ) continue;
					Real3  reldir = divide(relpos,dist);
					Real    vol_j = pData[setName_j]->vol[j];
					Real3    gWij = gW_ij(dist, reldir, smoothingLength);
					

          			Real3x3 v_j_gWij_rij = mult(- vol_j,tensorProduct(gWij,relpos));

					// Compute First derivative Renormalization Tensor
					L_i = add(L_i, v_j_gWij_rij);

					// Compute the reference Map Gradient
					pData[setName_i]->defoGrad[i][0] = add(pData[setName_i]->defoGrad[i][0], mult( (pData[setName_j]->originPos[j][0] - pData[setName_i]->originPos[i][0]) * vol_j, gWij));
					pData[setName_i]->defoGrad[i][1] = add(pData[setName_i]->defoGrad[i][1], mult( (pData[setName_j]->originPos[j][1] - pData[setName_i]->originPos[i][1]) * vol_j, gWij)); 
					pData[setName_i]->defoGrad[i][2] = add(pData[setName_i]->defoGrad[i][2], mult( (pData[setName_j]->originPos[j][2] - pData[setName_i]->originPos[i][2]) * vol_j, gWij));

				}
			}

			// Compute the inverse of the renormalization tensor.
			setDims(L_i,dims);
			L_i = inv3x3(L_i, -1);

			pData[setName_i]->defoGrad[i][0] = mult(pData[setName_i]->L[i],pData[setName_i]->defoGrad[i][0]);
			pData[setName_i]->defoGrad[i][1] = mult(pData[setName_i]->L[i],pData[setName_i]->defoGrad[i][1]);
			pData[setName_i]->defoGrad[i][2] = mult(pData[setName_i]->L[i],pData[setName_i]->defoGrad[i][2]);
 
			// // Invert the gradient of the reference map
			setDims(pData[setName_i]->defoGrad[i], dims);
			pData[setName_i]->defoGrad[i] = inv3x3(pData[setName_i]->defoGrad[i],-1);
			pData[setName_i]->defoGrad_withoutThermal[i] = mult(1.0 / thermalDeform, pData[setName_i]->defoGrad[i]);
			setDims(pData[setName_i]->defoGrad_withoutThermal[i], dims);

			// Compute the creep component			
			MatrixXd U = MatrixXd::Zero(3,3), R = MatrixXd::Zero(3,3);
			MatrixXd F_withoutThermal = MatrixXd::Zero(3,3); toMatrix3d(pData[setName_i]->defoGrad_withoutThermal[i], F_withoutThermal);

			polarDecompose(F_withoutThermal, R, U);						
			getStretch(U, pData[setName_i]->stretch_total[i], pData[setName_i]->stretch_dirs[i], pData[setName_i]->stretch_dirs_before[i]);



			toMatrix3d(pData[setName_i]->defoGrad_withoutThermal[i], U);
			pData[setName_i]->defoGrad_elastic[i] = pData[setName_i]->defoGrad_withoutThermal[i];

			pData[setName_i]->vol[i] = pData[setName_i]->originVol[i] * det3x3(pData[setName_i]->defoGrad[i]);
		}
	}
}


void SPHSolver::fixedPointIteration(Uint t){
	using namespace RealOps;

	currentTime += dt;

	std::cout << "   |----------------- Saving x(t)" << std::endl;
	for (const auto& setName : setNames){
		const int setID = ids[setName];
		const auto& ps = nsearch->point_set(setID);

		if(setName == "fluid"){
			#pragma omp parallel for num_threads(NUMTHREADS) 
			for (int i = 0; i < ps.n_points(); ++i){

				pData[setName]->posBefore[i]  = pData[setName]->pos[i];
				pData[setName]->velBefore[i]  = pData[setName]->vel[i];			
				pData[setName]->accBefore[i]  = pData[setName]->acc[i];			

				if(!pData[setName]->isThermalDirichlet[i]){
					pData[setName]->tempBefore[i]    = pData[setName]->temp[i];			
					pData[setName]->tempDotBefore[i] = pData[setName]->tempDot[i];			 
				} 

			}
		}
	}

	std::cout << "	|--- Placing Ghost Particles" << std::endl;
	addGhostParticles();
	reloadAndNeighborSearch();
	trimGhostParticles();

	for(Uint n = 0; n < (Uint) simData["fixedPointIterations"]; n ++){
		std::cout << "-------------------- Performing Fixed Point Iteration ... iteration " << n << std::endl;
		Real norm = 0;
		// computeDeformationGradient(t);
		computeDeformationGradient2(t);
		setGhostParticleTemperatures();
		computeInteractions(t);
		
		for (const auto& setName : setNames){

			const int setID = ids[setName];
			const auto& ps = nsearch->point_set(setID);

			if( setName == "fluid" ){
				#pragma omp parallel for num_threads(NUMTHREADS) 
				for (int i = 0; i < ps.n_points(); ++i){

					pData[setName]->pos[i] = add(
												pData[setName]->posBefore[i],
												 add(mult(dt, pData[setName]->velBefore[i]), 
													 mult(dt * dt / 4.0,add(pData[setName]->acc[i], pData[setName]->accBefore[i]))
													 )
												);	
					pData[setName]->pos[i][1] = pData[setName]->posBefore[i][1]; // Force zero deformation in the y-dir

					pData[setName]->vel[i] = add(
												pData[setName]->velBefore[i],
												mult(dt / 2.0, add(pData[setName]->acc[i], pData[setName]->accBefore[i]))
											);
					pData[setName]->vel[i][1] = 0; // Force zero deformation in the y-dir.

					if(!pData[setName]->isThermalDirichlet[i]){
						pData[setName]->temp[i] = pData[setName]->tempBefore[i] + (dt / 2.0) * (pData[setName]->tempDot[i] + pData[setName]->tempDotBefore[i]);
					}

					norm += length(pData["fluid"]->psi[i],pData[setName]->pos[i]);
					pData[setName]->psi[i] = pData[setName]->pos[i];

				}
			} else{
				// Do stuff for the boundary particles??
			}

		}
		
		std::cout << "dumb norm : " << norm << std::endl;
		// if(norm < 1.0E-9){
			// std::cout << "converged?" << std::endl;
			// break;
		// }
	}

	// XSPH(t);


	// for (const auto& setName : setNames){

	// 	const int setID = ids[setName];
	// 	const auto& ps = nsearch->point_set(setID);

	// 	if( setName == "fluid"){
	// 		#pragma omp parallel for num_threads(NUMTHREADS) 
	// 		for (int i = 0; i < ps.n_points(); ++i){

	// 			pData[setName]->vel[i] = add(
	// 									pData[setName]->vel[i],
	// 									mult(dt / 2.0, add(pData[setName]->acc[i], pData[setName]->accBefore[i]))
	// 									);

	// 		}
	// 	}
	// }

	
	
}




void SPHSolver::marchTime(Uint t){
	using namespace RealOps;

	std::cout << "-------------------- Timestep #" << t << std::endl;

	fixedPointIteration(t);

	// std::cout << "	|--- Updating Position " << std::endl;		

	// for (const auto& setName : setNames){

	// 	const int setID = ids[setName];
	// 	const auto& ps = nsearch->point_set(setID);


	// 	if( simData["computeElasticity"] == "Yes" ) computeDeformationGradient(t);
	// 	std::cout << "	|--- Calculating Interatctions" << std::endl;
	// 	computeInteractions(t);
	// 	std::cout << "	|--- Updating Vectors" << std::endl;

		

	// 	if( setName == "fluid"){
	// 	// Fluid Particles    : march the position, velocity, density.
	// 		#pragma omp parallel for num_threads(NUMTHREADS) 
	// 		for (int i = 0; i < ps.n_points(); ++i){
	// 			pData[setName]->vel[i]  = add(pData[setName]->vel[i], mult(dt,pData[setName]->acc[i]));
	// 			pData[setName]->pos[i]  = add(pData[setName]->pos[i], mult(dt,pData[setName]->vel[i]));

	// 			// ######################################
	// 			// Test For Deformation Gradient
	// 			// pData[setName]->pos[i] = add(pData[setName]-> pos[i] , Real3{(0.0 * pData[setName]->originPos[i][0] + 0.0  * pData[setName]->originPos[i][1]),
	// 																		//  (dt * 0.05 * pData[setName]->originPos[i][0] + 0.0  * pData[setName]->originPos[i][1]),0} );
	// 			// pData[setName]->pos[i] = add(pData[setName]-> pos[i] , Real3{dt * (0.5 * pData[setName]->originPos[i][0] + 0.0  * pData[setName]->originPos[i][1]),
	// 																		//  dt * (0.5 * pData[setName]->originPos[i][0] + 0.25 * pData[setName]->originPos[i][1]),0} );
	// 			// pData[setName]->pos[i] = Real3{(std::cos(currentTime) * pData[setName]->originPos[i][0] - std::sin(currentTime) * pData[setName]->originPos[i][1]),
	// 										//    (std::sin(currentTime) * pData[setName]->originPos[i][0] + std::cos(currentTime) * pData[setName]->originPos[i][1]),0};
	// 			// ######################################

	// 			// Density updates for WCSPH should only be considered for non-solid particles
	// 			if(!pData[setName]->isSolid[i]) pData[setName]->dens[i] = pData[setName]->dens[i] + dt * pData[setName]->densdot[i];
	// 			pData[setName]->temp[i] = pData[setName]->temp[i] + dt * pData[setName]->enthalpydot[i];					
	// 			pData[setName]->tau[i] = add(pData[setName]->tau[i],mult(dt, pData[setName]->tauDot[i]));				
	// 		}
	// 	} else if (setName == "boundary"){
	// 	// Boundary Particles : march the density only.
	// 		#pragma omp parallel for num_threads(NUMTHREADS)
	// 		for (int i = 0; i < ps.n_points(); ++i){
	// 			pData[setName]->dens[i] = pData[setName]->dens[i]  + dt * pData[setName]->densdot[i];
	// 			pData[setName]->temp[i] = pData[setName]->temp[i] + dt * pData[setName]->enthalpydot[i] / (pData[setName]->getSpecificHeat());
	// 		}
	// 	}
	// }


}


		// 	Unused code, formulation from Fatehi et al. Very sensitive with numerical precision.
		// 	Tensor4 Q_i, R_i;
		// 	Tensor3 S_i, P_i;
		// 	for (const auto& setName_j : setNames){

		// 		const int setID_j = ids[setName_j];
		// 		const auto& ps_j = nsearch->point_set(setID_j);

		// 		for (int _j = 0; _j < ps_i.n_neighbors(setID_j,i); _j++){

		// 			// Index of neighbor particle in setID_j
		// 			Uint const j = ps_i.neighbor(setID_j, i, _j);

		// 			// Define all the ingedients for the particle interaction
		// 			Real    rho_j = pData[setName_j]->dens[j];
		// 			Real3  relpos = sub(pData[setName_i]->pos[i],pData[setName_j]->pos[j]);
		// 			Real     dist = length(relpos);

		// 			if ( dist > smoothingLength) continue;
		// 			Real3  reldir = divide(relpos,dist);
		// 			Real    vol_j = pData[setName_j]->vol[j];
        //   			Real      Wij = W_ij(dist, smoothingLength);
		// 			Real3    gWij = gW_ij(dist, reldir, smoothingLength);

		// 			for (int m=0;m<3;m++) for(int n=0;n<3;n++) for(int o=0;o<3;o++) for(int p=0;p<3;p++)
		// 				R_i(m,n,o,p) = R_i(m,n,o,p) + vol_j * relpos[m] * reldir[n] * reldir[o] * gWij[p];

		// 			for (int m=0;m<3;m++) for(int n=0;n<3;n++) for(int k=0;k<3;k++)
		// 				S_i(m,n,k) = S_i(m,n,k) + vol_j * reldir[m] * reldir[n] * gWij[k];
					
		// 			for (int l=0;l<3;l++) for(int o=0;o<3;o++) for(int p=0;p<3;p++)
		// 				P_i(l,o,p) = P_i(l,o,p) + vol_j * relpos[l] * relpos[o] * gWij[p];
		// 		}
		// 	}


		// 	contract_323_to_4(S_i,L_i,P_i,Q_i);


		// 	Q_i.add(R_i);


			

		// 	for(int I=0;I<6;I++){
		// 		Uint m = IDXPAIR[I][0]; Uint n = IDXPAIR[I][1];
		// 		for(int J=0;J<6;J++){
		// 			Uint o = IDXPAIR[J][0]; Uint p = IDXPAIR[J][1];
		// 			_Q_i(I,J) = Q_i(m,n,o,p);
		// 		}
		// 	}
		// 	if(i == 1300){
		// 		std::cout << "before" << std::endl;
		// 		std::cout << _Q_i << std::endl;
		// 	}

		// 	if(dims == 2){
		// 		_Q_i(2,2) = 1.0; _Q_i(4,4) = 1.0; _Q_i(5,5) = 1.0;
		// 	} else if (dims == 1){
		// 		_Q_i(1,1) = 1.0; _Q_i(2,2) = 1.0; _Q_i(3,3) = 1.0; _Q_i(4,4) = 1.0; _Q_i(5,5) = 1.0;		
		// 	}

		// 	JacobiSVD<MatrixXd> svd_L2_i(_Q_i);
		// 	Real cond_second_i = svd_L2_i.singularValues()(0) / svd_L2_i.singularValues()(svd_L2_i.singularValues().size()-1);		
		// 	pData[setName_i]->isFS[i] = cond_second_i;			

		// 	L2_i = toReal3x3(_Q_i.fullPivLu().solve(neg_delta_mn));
		