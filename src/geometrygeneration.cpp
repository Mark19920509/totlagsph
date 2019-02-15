#include "geometrygeneration.hpp"

namespace GeometryGeneration{

    void randomSphere(Real3 center, Real r, Uint n, Real dx,
                      std::vector<Real>& x,
                      std::vector<Real>& y,
                      std::vector<Real>& z ){

		Real3 x_min = sub(center,mult(r,Real3{1.0,1.0,1.0}));
		Real3 x_max = add(center,mult(r,Real3{1.0,1.0,1.0}));
		std::cout << "... Generating Initial Configuration (Random Sphere) ..." << std::endl;

		std::vector<Real3> samples = thinks::poissonDiskSampling(0.8 * dx, x_min, x_max);

		int totalParticles = 0;
		for (int i = 0; i < samples.size(); i++ ){
			if(length(sub(samples[i],center)) > r) continue;
			x.push_back( samples[i][0] );
			y.push_back( samples[i][1] );
			z.push_back( samples[i][2] );
            totalParticles ++;
            if(totalParticles >= n) break;
		}
    }

    void randomLine(Real3 p1, Real3 p2, Uint n,
                    std::vector<Real>& x, 
                    std::vector<Real>& y, 
                    std::vector<Real>& z){
        
        Real l = p2[0] - p1[0];
        Real dx = l / (Real)n;
        Real smallestDist = 0.73 * dx;
        std::cout << "minimum stencil : " << smallestDist << std::endl;

        x.push_back(p1[0]);
        y.push_back(p1[1]);
        z.push_back(p1[2]);

        x.push_back(p2[0]);
        y.push_back(p2[1]);
        z.push_back(p2[2]);

        for (;;){
            std::mt19937 rng;
            rng.seed(std::random_device()());
            std::uniform_int_distribution<std::mt19937::result_type> randnum(0,100000);

            bool rejected = false;
            Real3 posToAdd = add(p1,Real3{ l *  ((Real)randnum(rng))/100000.0 ,0,0});

            for(int j=0;j<x.size();j++){
                Real3 neighPos = Real3{x[j],y[j],z[j]};                
                if ( length(sub(posToAdd,neighPos)) < smallestDist ){
                    rejected = true;
                    std::cout << "Rejected : " << length(sub(posToAdd,neighPos)) << std::endl;
                    break;
                }
            }
            if (rejected){
                continue;
            } else{
                x.push_back(posToAdd[0]);
                y.push_back(posToAdd[1]);
                z.push_back(posToAdd[2]);
                if (x.size() > n) break;
            }
        }       

    }

    void randomBox(Real3 p1, Real3 p2, Uint n, Real dx, 
                   std::vector<Real>& x, 
                   std::vector<Real>& y, 
                   std::vector<Real>& z){

		std::cout << "... Generating Initial Configuration (Random Box) ..." << std::endl;
        
		std::vector<Real3> samples = thinks::poissonDiskSampling(dx, p1, p2);
		int totalParticles = 0;
		for (;;){
            std::mt19937 rng;
            rng.seed(std::random_device()());
            std::uniform_int_distribution<std::mt19937::result_type> dist(0,samples.size()-1);
            
			x.push_back( samples[dist(rng)][0] );
			y.push_back( samples[dist(rng)][1] );
			z.push_back( samples[dist(rng)][2] );
            totalParticles ++;
            if(totalParticles >= n) break;
		}
        std::cout << "... Generation complete.." << std::endl;
    }

    void uniformCylinder(Real3 origin, Real length, Real radius, Uint nr, Uint nd, Real dx,
                        std::vector<Real>& x,
                        std::vector<Real>& y,
                        std::vector<Real>& z){
        Real pi = 3.141592;
        std::vector<Real3> pos;
        Real circum_thick = radius / (Real)nr;
        pos.push_back(origin);
        for (Uint i=1;i<=nr;i++){
            Uint ntheta = (Uint) (( 2.0 * pi * i * circum_thick) / circum_thick + 0.5);
            for (Uint j=0;j<ntheta;j++){
                Real dTheta = 2.0 * pi / ntheta;
                Real3 direction = Real3{0,cos(dTheta * j),sin(dTheta * j)};
                Real radius = (i * circum_thick);
                pos.push_back(add(origin,mult(radius,direction)));
            }
        }
        Real thick = (Real)(length / (Real)(nd-1));
        for (Uint layer = 0; layer < nd; layer ++){
            for (Uint i=0;i<pos.size();i++){
                x.push_back(pos[i][0] + thick * layer);
                y.push_back(pos[i][1]);
                z.push_back(pos[i][2]);
            }
        }
    }

    void randomCylinder(Real3 origin, Real length, Real radius, Uint n, Real dx,
                        std::vector<Real>& x, 
                        std::vector<Real>& y, 
                        std::vector<Real>& z){
        // Only a cylinder that shares the axis with the y axis works..
        Real3 x_min = sub(origin,mult(radius,Real3{0.0,1.0,1.0}));
        Real3 x_max = add(x_min,Real3{length, 2.0 * radius, 2.0 * radius});
        std::map<Uint,bool> didUse;
        std::cout << "... Generating Initial Configuration (Random Cylinder) ..." << std::endl;
                            
		std::vector<Real3> samples = thinks::poissonDiskSampling(0.8 * dx, x_min, x_max);
        int totalParticles = 0;

		for (int i = 0; i < samples.size(); i++ ){

            std::mt19937 rng;
            rng.seed(std::random_device()());
            std::uniform_int_distribution<std::mt19937::result_type> dist(0,samples.size()-1);
            Uint sampleIdx = dist(rng);
            Real3 sample = samples[sampleIdx];

            if ( sample[1] * sample[1] + sample[2] * sample[2] > radius * radius ) continue;
            if (didUse.count(sampleIdx) != 0) continue;

            didUse[sampleIdx] = true;

			x.push_back( sample[0] );
			y.push_back( sample[1] );
            z.push_back( sample[2] );

            totalParticles ++;
            if(totalParticles >= n) break;
		}
        
        for (int j = 0; j < x.size(); j++){

            if(x[j] < radius * 0.2){
                x[j] = 0.0;
            }            
            if(x[j] > length - radius * 0.2){
                x[j] = length;
            }
        }
        
            
    }

    void generateGeometry(json& simDataInput){

        std::vector<Real> x,y,z;
        std::cout << "... Generating Geometry." << std::endl;

        if(simDataInput["geometry"]["type"] == "sphere"){

            std::cout << "... Generating Random Sphere." << std::endl;
            Real dx = (Real) simDataInput["dx"];
            Uint n = (Uint) simDataInput["geometry"]["numParticles"];        
            Real3 center = Real3{(Real) simDataInput["geometry"]["dropletCenter"][0], (Real) simDataInput["geometry"]["dropletCenter"][1], (Real) simDataInput["geometry"]["dropletCenter"][2]};
            Real r = (Real) simDataInput["geometry"]["dropletRadius"];

            randomSphere(center,r,n,dx,x,y,z);

        } else if (simDataInput["geometry"]["type"] == "line"){

            std::cout << "... Generating Random Line." << std::endl;
            Real3 start = Real3{(Real) simDataInput["geometry"]["lineStart"][0],
                               (Real) simDataInput["geometry"]["lineStart"][1],
                               (Real) simDataInput["geometry"]["lineStart"][2]};

            Real3 end = Real3{(Real) simDataInput["geometry"]["lineEnd"][0],
                              (Real) simDataInput["geometry"]["lineEnd"][1],
                              (Real) simDataInput["geometry"]["lineEnd"][2]};

            Real dx = (Real) simDataInput["dx"];
            Uint n = (Uint) simDataInput["geometry"]["numParticles"];        
                      
            randomLine(start,end,n,x,y,z);

        } else if (simDataInput["geometry"]["type"] == "box"){

            std::cout << "... Generating Random Box." << std::endl;
            Real3 p1 = Real3{(Real)simDataInput["geometry"]["start"][0],(Real)simDataInput["geometry"]["start"][1],(Real)simDataInput["geometry"]["start"][2]};
            Real3 p2 = Real3{(Real)simDataInput["geometry"]["end"][0],(Real)simDataInput["geometry"]["end"][1],(Real)simDataInput["geometry"]["end"][2]};
            Real dx = (Real) simDataInput["dx"];
            Uint n = (Uint) simDataInput["geometry"]["numParticles"];        
    
            randomBox(p1,p2,n,dx,x,y,z);

            if (simDataInput["dimensions"] == 1){
                for (int i = 0; i < x.size(); i++){
                    y[i] = 0;
                    z[i] = 0;
                }
                x.push_back(0); y.push_back(0); z.push_back(0);
                x.push_back(p2[0]); y.push_back(0); z.push_back(0);
                
            } else if (simDataInput["dimensions"] == 2){
                for (int i = 0; i < x.size(); i++){
                    z[i] = 0;
                }
            }

        } else if (simDataInput["geometry"]["type"] == "uniform_cylinder"){
            std::cout << "... Generating Uniform Cylinder." << std::endl;            
            Real3 origin = Real3{(Real) simDataInput["geometry"]["origin"][0], (Real) simDataInput["geometry"]["origin"][1], (Real) simDataInput["geometry"]["origin"][2]};
            Real length = (Real) simDataInput["geometry"]["length"];
            Real radius = (Real) simDataInput["geometry"]["radius"];
            Uint n_rad = (Uint) simDataInput["geometry"]["nr"];
            Uint n_dir = (Uint) simDataInput["geometry"]["nd"];
            Real dx = (Real) simDataInput["dx"];
            
            uniformCylinder(origin, length, radius, n_rad, n_dir, dx, x, y, z);

        } else if (simDataInput["geometry"]["type"] == "cylinder"){

            std::cout << "... Generating Random Cylinder." << std::endl;            
            Real3 origin = Real3{(Real) simDataInput["geometry"]["origin"][0], (Real) simDataInput["geometry"]["origin"][1], (Real) simDataInput["geometry"]["origin"][2]};
            Real length = (Real) simDataInput["geometry"]["length"];
            Real radius = (Real) simDataInput["geometry"]["radius"];
            Real dx = (Real) simDataInput["dx"];
            Uint n = (Uint) simDataInput["geometry"]["numParticles"];        
    
            randomCylinder(origin, length, radius, n, dx, x, y, z);

        }

        std::string fname = simDataInput["geometry"]["geomOutputFile"].dump();
        fname = fname.substr(1, fname.size()); fname = fname.substr(0, fname.size()-1);
        
		H5PartFile* fileWriter = H5PartOpenFile(fname.c_str(),H5PART_WRITE);

        std::vector<double> doublex,doubley,doublez;
        for (int i=0;i<x.size();i++){
            doublex.push_back((double)x[i]);
            doubley.push_back((double)y[i]);
            doublez.push_back((double)z[i]);
        }
		H5PartSetStep(fileWriter,0);
		H5PartSetNumParticles(fileWriter, x.size());
		H5PartWriteDataFloat64(fileWriter,"x",&doublex[0]);
		H5PartWriteDataFloat64(fileWriter,"y",&doubley[0]);
		H5PartWriteDataFloat64(fileWriter,"z",&doublez[0]);
		H5PartCloseFile(fileWriter);
		std::cout << ">>> Geometry Generation Complete. Rejection Sampled to " << x.size() << " samples." << std::endl;


    }
    
}