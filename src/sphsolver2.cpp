#include "sphsolver.hpp"
#define CALC_HEAT 1

void SPHSolver::reloadAndNeighborSearch(){
	delete nsearch;
	nsearch = new NeighborhoodSearch((double)simData["smoothingLength"]*1.1,true);

	std::cout << "Realoading NSearch Algorithm" << std::endl;
	currentTime = 0.0;

	for (const auto& pDatEntry : pData){
		ids[pDatEntry.first] = nsearch->add_point_set( pDatEntry.second->pos.front().data(), pDatEntry.second->pos.size(), true, true);
		std::cout << "... Particle set \"" << pDatEntry.first << "\" with " << pDatEntry.second->numParticles << " Particles Set Loaded onto CompactNSearch." << std::endl;
		if (pDatEntry.second->pos.size() != pDatEntry.second->numParticles)	assert(!(pDatEntry.second->pos.size() == pDatEntry.second->numParticles));
	}


	totParticles = 0;
	Uint numSets = 0;

	for (const auto& pSet : nsearch->point_sets()){
		totParticles += pSet.n_points();
		numSets += 1;
	}

	std::cout << "... Loaded " << totParticles << " particles, with " << numSets << " point sets." << std::endl;
	std::cout << "... Searching Neighbors" << std::endl;

	auto t0 = std::chrono::high_resolution_clock::now();
	nsearch->find_neighbors();
	std::cout << "--- Neighborhood search took " << std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - t0).count() << "ms" << std::endl;

}

void SPHSolver::addGhostParticles(){
	using namespace RealOps;
	// const Real dx = (Real) simData["dx"];
	// const Uint dims = simData["dimensions"];
	// const Real smoothingLength = (Real) simData["smoothingLength"];
    const int m = (int) ((smoothingLength + (dx / 2.0)) / dx);


	json ghostDataInput = pData["ghost"]->parDataIn;
	delete pData["ghost"];
	pData["ghost"] = new ParticleAttributes(simData, ghostDataInput);
	
    ghostMap.clear();

	for (const auto& setName_i : setNames){

		const int setID_i = ids[setName_i];
		const auto& ps_i = nsearch->point_set(setID_i);
        if (setName_i != "fluid") continue;

		for (int i = 0; i < ps_i.n_points(); ++i){

			if ( pData[setName_i]->particleDensity[i] > simData["ghostTol"]) continue;		            
            // std::cout << pData[setName_i]->particleDensity[i] << std::endl;
			Real3& pos_i = pData[setName_i]->pos[i];
            for (int ii = -m; ii <= m; ii ++){
                
            for (int jj = -m; jj <= m; jj ++){
                if (dims == 2){

                    int nx = ((int)(pos_i[0] / dx) - ii);
                    int ny = ((int)(pos_i[1] / dx) - jj);
                    int nz = 0;
                    std::tuple<int,int,int> key = std::make_tuple(nx,ny,nz);

                    // Ghost particle is already filled
                    if (ghostMap.count(key) > 0) continue;                    
                    ghostMap[key] = true;
					// std::cout << "jarjar" << std::endl;
                    Real3 posToPush{nx * dx, ny * dx, 0};                        
                    pData["ghost"]->addDefaultFluidParticleAtPosition(posToPush);
					// std::cout << "binks" << std::endl;
                    pData["ghost"]->isGhost.back() = true;
					// std::cout << "binks" << std::endl;
                }
            }}

            
		}
	}    

}

void SPHSolver::trimGhostParticles(){

	using namespace RealOps;
	// const Real dx = (Real) simData["dx"];
	// const Uint dims = simData["dimensions"];
	// const Real smoothingLength = (Real) simData["smoothingLength"];
    const int m = (int) ((smoothingLength + (dx / 2.0)) / dx);

	for (const auto& setName_i : setNames){

		const int setID_i = ids[setName_i];
		const auto& ps_i = nsearch->point_set(setID_i);
		std::cout << setName_i << std::endl;
        if (setName_i != "ghost") continue;

		#pragma omp parallel for num_threads(NUMTHREADS) 

		for (int i = 0; i < ps_i.n_points(); ++i){

			pData[setName_i]->particleDensity[i] = 0;

			for (const auto& setName_j : setNames){

				const int setID_j = ids[setName_j];
				const auto& ps_j = nsearch->point_set(setID_j);

				for (int _j = 0; _j < ps_i.n_neighbors(setID_j,i); _j++){

					Uint const j = ps_i.neighbor(setID_j, i, _j);
					Real3  relpos = sub(pData[setName_i]->pos[i],pData[setName_j]->pos[j]);
					Real     dist = length(relpos);	
					Real      Wij = W_ij(dist, smoothingLength);					

					if (dist < smoothingLength){
						pData[setName_i]->particleDensity[i] += pData[setName_j]->vol[j] * Wij; 
					}

				}

			}	
			pData[setName_i]->particleDensity[i] += pData[setName_i]->vol[i] * W_ij(0, smoothingLength);
			if (pData[setName_i]->particleDensity[i] > (Real)simData["ghostTrimTol"]){
				pData[setName_i]->isActive[i] = false;
			}
		}
	}    	
}

void SPHSolver::setGhostParticleTemperatures(){
	std::cout << "		|------ ... Setting Ghost Particle Temperatures"<< std::endl;

	for (const auto& setName_i : setNames){

		const int setID_i = ids[setName_i];
		const auto& ps_i = nsearch->point_set(setID_i);

		#pragma omp parallel for num_threads(NUMTHREADS)
		for (int i = 0; i < ps_i.n_points(); ++i){
			bool isBoundary_i = (setName_i == "boundary") ? true : false;
			bool isGhost_i = (setName_i == "ghost") ? true : false;

			// Only need to set temperature for the ghost particles
			if(!isGhost_i) continue;

			for (const auto& setName_j : setNames){
				const int setID_j = ids[setName_j];
				const auto& ps_j = nsearch->point_set(setID_j);

				const bool isBoundary_j = (setName_j == "boundary") ? true : false;
				const bool isGhost_j = (setName_j == "ghost") ? true : false;				

				// If the neighbor is another ghost particle, we don't care.
				if(isGhost_j) continue;
				// Loop over the neighbors of particl i in the indexed particle set.
				Real minDist = 10 * dx;
				Real minDistTemp = -1;

				for (int _j = 0; _j < ps_i.n_neighbors(setID_j,i); _j++){
					
					// Index of neighbor particle in setID_j
					Uint const j = ps_i.neighbor(setID_j, i, _j);
					// Define all the ingedients for the particle interaction
					Real3  relpos = sub(pData[setName_i]->pos[i],pData[setName_j]->pos[j]);
					Real3 relposj = mult(-1.0,relpos);

					Real     dist = length(relpos);
					if ( dist > smoothingLength) continue;

					if ( dist < minDist ){
						minDist = dist;
						minDistTemp = pData[setName_j]->temp[j];
					}

				}
				pData[setName_i]->temp[i] = minDistTemp;
			}
		}
	}
}


void SPHSolver::computeInteractions(Uint t){
	using namespace RealOps;

	// The fluid.
	// const Real dx = (Real) simData["dx"];
	// const Real dt = (Real) simData["dt"];
	// const Uint dims = simData["dimensions"];
	// const Real smoothingLength = (Real) simData["smoothingLength"];

	std::cout << "		|------ ... Performing Fluid - (Fluid / Boundary) Interactions"<< std::endl;
	for (const auto& setName_i : setNames){

		const bool isGhost_i = (setName_i == "ghost");
		const int setID_i = ids[setName_i];
		const auto& ps_i = nsearch->point_set(setID_i);

		// Ghost particles don't have any time derivatives to compute
		if (isGhost_i) continue;

		#pragma omp parallel for num_threads(NUMTHREADS)

		for (int i = 0; i < ps_i.n_points(); ++i){

			const Real&   m_i = pData[setName_i]->mass[i];
			const Real& rho_i = pData[setName_i]->dens[i];
			const Real&   T_i = pData[setName_i]->temp[i];
      		const Real& vol_i = pData[setName_i]->vol[i];
			const Real&   k_i = thermalConductivity(pData[setName_i]->getType(),
				  							 T_i);

			const Real& vol_i_o = pData[setName_i]->originVol[i];
			const Real&   P_i = EOS(T_i, pData[setName_i]->getT0(),
									     rho_i, pData[setName_i]->getRho0(),
									     pData[setName_i]->getSoundSpeed(),
										 pData[setName_i]->getThermalExpansion());

			const Real3& vel_i   = pData[setName_i]->vel[i];
			const Real3& pos_i_o = pData[setName_i]->originPos[i];
			const Real3x3& L_i  = pData[setName_i]->L[i];
			const Real& pDens_i = pData[setName_i]->particleDensity[i];

			bool isSensor_i   = (pData[setName_i]->isSensor[i]);
			bool isBoundary_i = (setName_i == "boundary") ? true : false;
			bool isGhost_i = (setName_i == "ghost") ? true : false;

			pData[setName_i]->tempGrad[i] = zerovec;
			pData[setName_i]->acc[i] = zerovec;
			pData[setName_i]->tempDot[i] = 0;


			for (const auto& setName_j : setNames){

				const int setID_j = ids[setName_j];
				const auto& ps_j = nsearch->point_set(setID_j);
				const bool isBoundary_j = (setName_j == "boundary") ? true : false;
				const bool isGhost_j = (setName_j == "ghost") ? true : false;				
				
				// Loop over the neighbors of particl i in the indexed particle set.
				for (int _j = 0; _j < ps_i.n_neighbors(setID_j,i); _j++){

					// Index of neighbor particle in setID_j
					Uint const j = ps_i.neighbor(setID_j, i, _j);

					// Skip if the particle is not active.
					// std::cout << j << std::endl;
					// const bool isGhostActive_j = (pData[setName_j]->isGhost[j] && pData[setName_j]->isActive[j]);
					// if(!isGhostActive_j) continue;
					if(!pData[setName_j]->isActive[j]) continue;
					// std::cout << "111" << std::endl;

					// Define all the ingedients for the particle interaction
					Real3  relpos = sub(pData[setName_i]->pos[i],pData[setName_j]->pos[j]);
					Real3 relposj = mult(-1.0,relpos);

					Real     dist = length(relpos);
					if ( dist > smoothingLength) continue;
					
					const Real&      m_j   = pData[setName_j]->mass[j];
					const Real&    vol_j   = pData[setName_j]->vol[j];

					const Real3& vel_j = pData[setName_j]->vel[j];
					const Real3& normal_j =  pData[setName_j]->normalVec[j];
					const Real& pDens_j = pData[setName_j]->particleDensity[j];

					Real   rho_i_temp, rho_j_temp;
					Real3  relvel = sub(pData[setName_i]->vel[i],pData[setName_j]->vel[j]);
					Real3  reldir = divide(relpos,dist);

					Real      Wij =  W_ij(dist,         smoothingLength);
					Real3    gWij = gW_ij(dist, reldir, smoothingLength);


					Real    rho_j = pData[setName_j]->dens[j];
					Real3 densGrad_j = pData[setName_j]->densGrad[j];
					Real   	  T_j = pData[setName_j]->temp[j];

					// std::cout << "wagya" << std::endl;					
					// std::cout << "wagya!" << std::endl;

					// Actual density of the boundary material, to account for the correct thermal diffusivity.
					rho_i_temp = isBoundary_i ? pData[setName_i]->getMaterialDensity() : rho_i ;
					rho_j_temp = isBoundary_j ? pData[setName_j]->getMaterialDensity() : rho_j ;

					// If the laplacian corrector cannot be defined, use the conventional operator.						
					Real heat;
					heat = inconsistentHeatTransfer(T_i, T_j, m_j, rho_i_temp, rho_j_temp,
													relpos, reldir,
													dist, gWij, vol_j);		

					pData[setName_i]->tempGrad[i] = add(pData[setName_i]->tempGrad[i], 
														mult((T_j-T_i)*vol_j, gWij));
					// pData[setName_i]->enthalpydot[i] += heat * 15.0 / (7900.0 * 450.0);
					pData[setName_i]->tempDot[i] += heat;

					// Damping for solid model
					if( !isGhost_i && !isGhost_j){

						pData[setName_i]->acc[i] = add(pData[setName_i]->acc[i],
														mult(	- pData[setName_j]->mass[j] / pData[setName_i]->dens[i],
																mult(
																	scalarViscosity_ij( (Real) simData["Cl"],
																						(Real) simData["Cq"],
																							smoothingLength,
																						(Real) pData[setName_i]->getSoundSpeed(),
																							relpos,
																							relvel
																					),
																	gWij
																)
														)
							 						   );

						pData[setName_i]->acc[i] = add(pData[setName_i]->acc[i],
								mult(	- (Real) simData["damping"] * pData[setName_j]->vol[j],
										mult(Wij, relvel))
								);

					}

				}

			}

			// If particle i is not a boundary particle, accumulate the momentum contribution.
			if (!isBoundary_i){
				// Loop over the initial neighbors.
				for(const auto& j_ : pData[setName_i]->nMap[i]){
					// If i is a solid particle and j is a solid particle, accumulate the cauchy stress contributions
					const std::string setName_j = std::get<0>(j_);
					const Uint j = std::get<1>(j_);

					if ( pData[setName_i]->isSolid[i] ){
						if( pData[setName_j]->isSolid[j] ){

							const Real&    vol_j_o = pData[setName_j]->originVol[j];
							const Real3&   pos_j_o = pData[setName_j]->originPos[j];
							Real3  relpos_o = sub(pos_i_o,pos_j_o);
							Real3  gWij_o = gW_ij(length(relpos_o), dir(relpos_o), smoothingLength);

							// from ganzenmuller
							pData[setName_i]->acc[i] = add(pData[setName_i]->acc[i], 			
															mult(   1.0 * vol_i_o * vol_j_o / m_i,
																	add(											
																		mult(mult(pData[setName_i]->defoGrad[i],pData[setName_i]->secondPKStress[i]), mult(pData[setName_i]->L_o[i], gWij_o)),
																		mult(mult(pData[setName_j]->defoGrad[j],pData[setName_j]->secondPKStress[j]), mult(pData[setName_j]->L_o[j], gWij_o))
																	)
																)
														  );

						}
						else{
							// Implement this part for FSI (non-solid fluid (j) <> solidified fluid (i) )
						}
					}

				}
			}


			if (isBoundary_i){
				pData[setName_i]->acc[i] = zerovec;
				pData[setName_i]->vel[i] = zerovec;
			} else{
				// pData[setName_i]->acc[i] = add(pData[setName_i]->acc[i],Real3{0.0,-0.2,0.0});
				//prescribed boundary
				// if(pData[setName_i]->originPos[i][0] >= 1.9 && pData[setName_i]->originPos[i][1] >= 0.474 ){
				// 	pData[setName_i]->acc[i][1] = (-10.0 * std::min((Real)t/(100.0),(Real)1.0));
				// 	// pData[setName_i]->pos[i][1] = pData[setName_i]->originPos[i][1] - std::min(0.1 * t / (200.0)  ,0.1);
				// }
			}
			
		}

	}
}





void SPHSolver::computeDeformationGradient2(Uint t){
	using namespace RealOps;



	for (const auto& setName_i : setNames){

		const int setID_i = ids[setName_i];
		if (setName_i != "fluid") continue;
		const auto& ps_i = nsearch->point_set(setID_i);
		
		#pragma omp parallel for num_threads(NUMTHREADS)
		for (int i = 0; i < ps_i.n_points(); ++i){
			Real3x3& L_i_o = pData[setName_i]->L_o[i];			
			Real3x3& F_i = pData[setName_i]->defoGrad[i];
			Real3x3& F_theta_i   = pData[setName_i]->defoGrad_thermal[i];
			Real3x3& F_elastic_i = pData[setName_i]->defoGrad_elastic[i]; 

			Real& vol_i = pData[setName_i]->vol[i];
			Real3x3& secondPKStress_i = pData[setName_i]->secondPKStress[i];
			secondPKStress_i = zeromat;

			const Real&    vol_i_o = pData[setName_i]->originVol[i];
			const Real3&   pos_i_o = pData[setName_i]->originPos[i];
			const Real3&   pos_i = pData[setName_i]->pos[i];
			const Real&   temp_i = pData[setName_i]->temp[i];

			F_i = zeromat;
			F_theta_i = add(identity(), mult(pData[setName_i]->getThermalExpansion() * (Real) temp_i, identity()));
			if(t == 0){ L_i_o = zeromat; }


			for(const auto& j_ : pData[setName_i]->nMap[i]){
				const std::string setName_j = std::get<0>(j_);
				const Uint j = std::get<1>(j_);
				if (setName_j != "fluid") continue;
				const Real&    vol_j_o = pData[setName_j]->originVol[j];
				const Real3&   pos_j_o = pData[setName_j]->originPos[j];
				const Real3&   pos_j = pData[setName_j]->pos[j];

				Real3  relpos_o = sub(pos_i_o,pos_j_o);
				Real     dist_o = length(relpos_o);	
				Real      Wij_o = W_ij(dist_o, smoothingLength);					
				Real3    gWij_o = gW_ij(length(relpos_o), dir(relpos_o), smoothingLength);

				if (dist_o > smoothingLength) continue;
		
				// Compute the renormalization tensor for the initial configuration.
				if(t == 0){
					Real3x3 v_j_gWij_rij = mult(- vol_j_o,tensorProduct(gWij_o,relpos_o));
					L_i_o = add(L_i_o, v_j_gWij_rij);
				}

				F_i[0] = add(F_i[0], mult( (pos_j[0]- pos_i[0]) * vol_j_o, gWij_o));
				F_i[1] = add(F_i[1], mult( (pos_j[1]- pos_i[1]) * vol_j_o, gWij_o)); 
				F_i[2] = add(F_i[2], mult( (pos_j[2]- pos_i[2]) * vol_j_o, gWij_o));
				
			}

			if (t == 0){ setDims(L_i_o,dims); L_i_o = inv3x3(L_i_o, -1); }
			F_i[0] = mult(L_i_o,F_i[0]); F_i[1] = mult(L_i_o,F_i[1]); F_i[2] = mult(L_i_o,F_i[2]); 

			setDims(F_i, dims);
			
			F_elastic_i = mult(F_i, inv3x3(F_theta_i, dims));
			// print(F_i);
			// print(F_theta_i);
			// print(F_elastic_i);


			if( setName_i == "fluid" ){	
				// Infinitesimal strain
				Real3x3 strainTensor_i = mult(0.5, sub( mult(transpose(F_elastic_i),F_elastic_i), identity()));
				Real3x3 strainTensorVol_i = mult(0.333333333333 * trace(strainTensor_i), identity());
				Real3x3 strainTensorDev_i = sub(strainTensor_i, strainTensorVol_i);

				// Second Piola-Kirchhoff Stress				
				secondPKStress_i = add(
										mult(     (Real) pData[setName_i]->getFirstLame() , strainTensorVol_i),
										mult(2. * (Real) pData[setName_i]->getShearModulus(), strainTensorDev_i)
									);
			}

			vol_i = vol_i_o * det3x3(F_i);
		}
	}
}