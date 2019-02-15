

class ParticleData{

	int numParticles;
	std::vector<std::array<Real,3>> pos();
	std::vector<std::array<Real,3>> vel();
	std::vector<std::array<Real,3>> acc();
	std::vector<std::array<Real>> dens();
	std::vector<std::array<Real>> temp();

	ParticleData(){

	}


};