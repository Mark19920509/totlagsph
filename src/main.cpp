#include <CompactNSearch>

#include <iostream>
#include <fstream>
#include "json.hpp"
#include "helpers.hpp"
#include "sphsolver.hpp"
#include "particleAttrib.hpp"
#include "filewriter.hpp"
#include "geometrygeneration.hpp"

using json = nlohmann::json;
using namespace CompactNSearch;

int main ( void ){

	std::cout << ">>> Starting Program." << std::endl;
	std::map<std::string,ParticleAttributes*> pData;
	json simDataInput, fluidDataInput, boundaryDataInput, ghostDataInput;
	std::cout << ">>> Reading system parameters." << std::endl;
	readJSON(simDataInput,      "../inputs/constants.json");
	std::cout << ">>> Reading fluid Data." << std::endl;
	readJSON(fluidDataInput,    "../inputs/fluid.json");
	std::cout << ">>> Reading boundary Data." << std::endl;
	readJSON(boundaryDataInput, "../inputs/boundary.json");
	std::cout << ">>> Reading Ghost Particle Data." << std::endl;
	readJSON(ghostDataInput, "../inputs/ghost.json");

	if( simDataInput["operation"] == "SPHSimulation" ){
		
		std::cout << ">>> Operation : SPHSimulation" << std::endl;

		pData["fluid"] = new ParticleAttributes(simDataInput, fluidDataInput);
		pData["ghost"] = new ParticleAttributes(simDataInput, ghostDataInput);

		// std::cout << ">>> End of Particle Definition." << std::endl;

		FileWriter writer(simDataInput["outputFile"].dump());

		SPHSolver solver(simDataInput, pData, 
						(Real) simDataInput["dx"],
						(Real) simDataInput["dt"],
						(Real) simDataInput["smoothingLength"],
						(Uint) simDataInput["dimensions"] );

		std::cout << ">>> Writing Initial Configuration (t=0)" << std::endl;
		writer.write(pData,simDataInput["smoothingLength"]);

		std::cout << ">>> Timestep start" << std::endl;

			for(int t=0;t<(Uint)simDataInput["steps"];t++){
				
				solver.marchTime(t);

				if (t >= (Uint)simDataInput["outputStart"]){
					// std::cout << "asdf : " << (t+1) % ((Uint)simDataInput["outper"]) << std::endl;
				if ((t+1) % ((Uint)simDataInput["outper"]) == 0)
					writer.write(pData,simDataInput["smoothingLength"]);
				if ((t+1) % (Uint)simDataInput["newfileper"] == 0)
					writer.openNextFile();
				}

				std::cout << "End of timestep : " << t << std::endl;
			}

		writer.close();

	} else if ( simDataInput["operation"] == "GeometryGeneration" ){
		std::cout << "MODE: Geometry Generation mode." << std::endl;
		GeometryGeneration::generateGeometry(simDataInput);

	}

	std::cout << ">>> Exiting Program." << std::endl;
	return 0;

}
