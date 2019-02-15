#ifndef __HELPERS__
#define __HELPERS__

#include <iostream>
#include <fstream>
#include <CompactNSearch>
#include <Eigen/Dense>
#include <math.h>
#include "json.hpp"
#include "vvector.h"

#define NUMTHREADS 4
#define MAX_NEIGHBORS 500
#define EPSL_SMALL 1.0E-10
#define EPSL_SMALL2 1.0E-15
#define PI 3.1415926535897932384626433832795028841971693993751058
// using namespace CompactNSearch;
using json = nlohmann::json;
// using Real = double;

inline void readJSON(json& jsonobj, std::string fname){
	std::ifstream input(fname);
	input >> jsonobj;
}

namespace RealOps{
	using Uint = unsigned int;
	using namespace CompactNSearch;
	using namespace Eigen;

	using Real4      = std::array<Real,4>;

	using Real3      = std::array<Real,3>;
	using Real3x3    = std::array<Real3,3>;
	using Real3x3x3  = std::array<Real3x3,3>;

	using Real6      = std::array<Real,6>;
	using Real6x6    = std::array<Real6,6>;
	// (0,0) (0,1) (0,2)
	// (1,0) (1,1) (1,2)
	// (2,0) (2,1) (2,2)
	
	// const Uint IDXPAIR_2D[4][2] = {{0,0},{1,1},{0,1}};
	const Uint IDXPAIR_2D[4][2] = {{0,0},{1,1},{0,1}};
	//                             
	const Uint IDXPAIR[6][2] = {{0,0},{1,1},{2,2},{0,1},{1,2},{0,2}};

	const Uint IDXPAIR_FULL_1D[3][3] = {{0,1,2},{3,4,5},{6,7,8}};
	// const Real NEG_DELTA_MN[6]   = {-1.0,-1.0,-1.0,0,0,0};
	// Eigen::Matrix3d A((Eigen::Matrix3d() << 1, 2, 3, 4, 5, 6, 7, 8, 9).finished());

	inline Real3x3 toReal3x3 (VectorXd vec){
		return Real3x3{Real3{vec(0),vec(3),vec(5)},
					   Real3{vec(3),vec(1),vec(4)},
					   Real3{vec(5),vec(4),vec(2)}};
	}

	inline void toMatrix3d (Real3x3& b,MatrixXd& a){
		
		for(int i=0;i<3;i++) for(int j=0;j<3;j++)
			a(i,j) = b[i][j];
	}

	inline void toReal3x3 (MatrixXd& b,Real3x3& a){
		
		for(int i=0;i<3;i++) for(int j=0;j<3;j++)
			a[i][j] = b(i,j);
		
	}

	inline void toReal3 (MatrixXd b, Real3& a){
		
		for(int i=0;i<3;i++)
			a[i] = (Real)b(i);
		
	}

	// neg_delta_mn(0) = -1.0;	neg_delta_mn(1) = -1.0;	neg_delta_mn(2) = -1.0;
	// neg_delta_mn(3) =  0.0; neg_delta_mn(4) =  0.0; neg_delta_mn(5) =  0.0;
	// IDXPAIR[6][2] = {{0,0},{1,1},{2,2},{0,1},{1,2},{0,2}};

	inline Real3x3 toReal3x3From6 (VectorXd& vec, Uint dims){

		if(dims == 2){
			vec(2) = 1.0;
		} else if (dims == 1){
			vec(1) = 1.0; vec(2) = 1.0;
		}

		return Real3x3{Real3{vec(0),vec(3),vec(5)},
					   Real3{vec(3),vec(1),vec(4)},
					   Real3{vec(5),vec(4),vec(2)}};
	}


	inline Real maxEig( VectorXd& a){
		Real result = a[0];
		if (a[1] > result)
		  result = a[1];

		if (a[2] > result)
		  result = a[2];

		return result;
	}

	inline Real minEig( VectorXd& a){
		Real result = a[0];
		if (a[1] < result)
		  result = a[1];

		if (a[2] < result)
		  result = a[2];

		return result;
	}

	inline void assign( MatrixXd& a, const Real3x3& b ){
		for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			a(i,j)=b[i][j];
	}

	inline void assign( Real3x3& b, const MatrixXd& a ){
		for(int i=0;i<3;i++)
		for(int j=0;j<3;j++)
			b[i][j] = a(i,j);
	}


	inline std::array<Real,3> add( const std::array<Real,3>& a,  const std::array<Real,3>& b) {
	    return std::array<Real,3>{a[0]+b[0],a[1]+b[1],a[2]+b[2]};
	}


	inline Real3 sub( Real3 a,  Real3 b) {
		return Real3{a[0]-b[0],a[1]-b[1],a[2]-b[2]};
	}

	inline Real3x3 sub( Real3x3 a,  Real3x3 b) {
		return Real3x3{sub(a[0],b[0]),
					   sub(a[1],b[1]),
					   sub(a[2],b[2])};
	}

	// inline Real3x3 sub( Real3x3 a,  Real3x3 b) {
	// 	return Real3x3{Real3{a[0][0] - b[0][0],a[0][1] - b[0][1],a[0][2] - b[0][2]},
	// 				   Real3{a[1][0] - b[1][0],a[1][1] - b[1][1],a[1][2] - b[1][2]},
	// 				   Real3{a[2][0] - b[2][0],a[2][1] - b[2][1],a[2][2] - b[2][2]}};
	// }

	inline Real3 mult( Real c, Real3 a) {
	    return Real3{c * a[0], c * a[1], c * a[2]};
	}

	inline Real3x3 mult( Real3x3x3 c, Real3 a) {
	    return Real3x3{Real3{c[0][0][0]*a[0] + c[0][0][1]*a[1] + c[0][0][2]*a[2], 
							 c[0][1][0]*a[0] + c[0][1][1]*a[1] + c[0][1][2]*a[2], 
							 c[0][2][0]*a[0] + c[0][2][1]*a[1] + c[0][2][2]*a[2]},

					   Real3{c[1][0][0]*a[0] + c[1][0][1]*a[1] + c[1][0][2]*a[2], 
							 c[1][1][0]*a[0] + c[1][1][1]*a[1] + c[1][1][2]*a[2], 
							 c[1][2][0]*a[0] + c[1][2][1]*a[1] + c[1][2][2]*a[2]},

					   Real3{c[2][0][0]*a[0] + c[2][0][1]*a[1] + c[2][0][2]*a[2], 
							 c[2][1][0]*a[0] + c[2][1][1]*a[1] + c[2][1][2]*a[2], 
							 c[2][2][0]*a[0] + c[2][2][1]*a[1] + c[2][2][2]*a[2]}};
	}

	inline Real3 mult( Real3x3x3 a, Real3x3 b) {
		Real3 c{0,0,0};
		for (int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++) c[i] = c[i] + a[i][j][k] * b[j][k];
	    return c;
	}


	// cij = aik*bkj
	inline Real3x3 mult( const Real3x3& a, const Real3x3& b) {
	    return Real3x3{Real3{a[0][0]*b[0][0]+a[0][1]*b[1][0]+a[0][2]*b[2][0], 
							 a[0][0]*b[0][1]+a[0][1]*b[1][1]+a[0][2]*b[2][1], 
							 a[0][0]*b[0][2]+a[0][1]*b[1][2]+a[0][2]*b[2][2]},

					   Real3{a[1][0]*b[0][0]+a[1][1]*b[1][0]+a[1][2]*b[2][0], 
					         a[1][0]*b[0][1]+a[1][1]*b[1][1]+a[1][2]*b[2][1], 
							 a[1][0]*b[0][2]+a[1][1]*b[1][2]+a[1][2]*b[2][2]},

					   Real3{a[2][0]*b[0][0]+a[2][1]*b[1][0]+a[2][2]*b[2][0], 
					   		 a[2][0]*b[0][1]+a[2][1]*b[1][1]+a[2][2]*b[2][1], 
							 a[2][0]*b[0][2]+a[2][1]*b[1][2]+a[2][2]*b[2][2]}};
	}

	// inline Real3x3 mult( Real c, const Real3x3& a) {
	//     return Real3x3{mult(c,a[0]),
	// 				   mult(c,a[1]),
	// 				   mult(c,a[2])};
	// }

	inline std::array<Real,3> divide(const std::array<Real,3>& a, Real c) {
	    return std::array<Real,3>{a[0] / c, a[1] / c, a[2] / c};
	}

	inline std::array<Real,3> neg(const std::array<Real,3>& a) {
	    return std::array<Real,3>{ -a[0], -a[1], -a[2]};
	}

	inline std::array<Real,3> dir( const std::array<Real,3>& a ) {
	    return divide(a,a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
	}

	inline std::array<Real,3> cross( const std::array<Real,3>& a,  const std::array<Real,3>& b) {
		return std::array<Real,3>{ a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]};
	}

	inline Real dot( const std::array<Real,3>& a,  const std::array<Real,3>& b) {
	    return (a[0] * b[0] + a[1] * b[1] + a[2] * b[2]);
	}

	inline Real3 normalize( const Real3& a ){
		Real l = std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]) + EPSL_SMALL;
		if ( l < 1.E-4){
			return Real3{0,0,0};
		} else{
		return Real3{ a[0]/l,
					  a[1]/l,
					  a[2]/l
					};
		}
	}

	inline Real length( const std::array<Real,3>& a ) {
	    return std::sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
	}

	inline Real length( const std::array<Real,3>& a , const std::array<Real,3>& b) {
	    return std::sqrt( (a[0]-b[0]) * (a[0]-b[0]) + (a[1]-b[1]) * (a[1]-b[1]) + (a[2]-b[2]) * (a[2]-b[2]));
	}

	inline Real length2( const std::array<Real,3>& a ) {
	    return (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
	}

	inline bool iseq( const std::array<Real,3>& a, const std::array<Real,3>& b){
		return (a[0]==b[0] && a[1]==b[1] && a[2]==b[2]) ? true : false;
	}
	inline void print( const Real3x3& a){
		std::cout << "\033[1;31m -------- Real3x3 -------- \033[0m" << std::endl;
		std::cout << " Real3 : { " << a[0][0] << ", " << a[0][1] << ", " << a[0][2] << " }" << std::endl;
		std::cout << " Real3 : { " << a[1][0] << ", " << a[1][1] << ", " << a[1][2] << " }" << std::endl;
		std::cout << " Real3 : { " << a[2][0] << ", " << a[2][1] << ", " << a[2][2] << " }" << std::endl;
		
	}

	inline void print(  const std::string& s , const std::array<Real,3>& a ){
		std::cout << s << " Real3 : { " << a[0] << ", " << a[1] << ", " << a[2] << " }" << std::endl;
 	}

	inline void print(  const std::string& s , const Real& a){
		std::cout << s << " Real : { " << a << " }"  << std::endl;
 	}


 	inline Real3x3 mult(Real a, Real3x3 b){
 		return Real3x3{mult(a,b[0]),mult(a,b[1]),mult(a,b[2])};
 	}

 	inline Real3 mult(Real3x3 a, Real3 b){
 		return Real3{ a[0][0] * b[0] + a[0][1] * b[1] + a[0][2] * b[2],
		 			  a[1][0] * b[0] + a[1][1] * b[1] + a[1][2] * b[2],
					  a[2][0] * b[0] + a[2][1] * b[1] + a[2][2] * b[2]};
 	}



	inline Real3x3 add( Real3x3 a,  Real3x3 b) {
	    return Real3x3 { add(a[0],b[0]), add(a[1],b[1]), add(a[2],b[2])};
	}

	inline Real3x3x3 add( Real3x3x3 a,  Real3x3x3 b) {
	    return Real3x3x3 { add(a[0],b[0]), add(a[1],b[1]), add(a[2],b[2])};
	}

	inline Real3x3 transpose( const Real3x3 a){
		// return Real3x3{ Real3{a[0][0],a[0][1],a[0][2]},
		// 				   Real3{a[1][0],a[1][1],a[1][2]},
		// 				   Real3{a[2][0],a[2][1],a[2][2]}};

		return Real3x3{ Real3{a[0][0],a[1][0],a[2][0]},
						Real3{a[0][1],a[1][1],a[2][1]},
						Real3{a[0][2],a[1][2],a[2][2]}};
	}

	inline Real doubleDot(const Real3x3& a, const Real3x3& b){
		Real bar=0.0;
		for(int i=0;i<3;i++)for(int j=0;j<3;j++)
			bar+=(a[i][j]*b[i][j]);
		return bar;
	}
	// inline Real3x3 add( const Real3x3& a,  const Real3x3& b) {
	//     return Real3x3 { Real3{a[0][0] + b[0][0],a[0][1] + b[0][1],a[0][2] + b[0][2]},
	// 				     Real3{a[1][0] + b[1][0],a[1][1] + b[1][1],a[1][2] + b[1][2]},
	// 				     Real3{a[2][0] + b[2][0],a[2][1] + b[2][1],a[2][2] + b[2][2]} };
	// }

	inline Real6 toVoigt(Real3x3& a, Uint dims){
		return Real6{a[0][0], a[1][1], a[2][2], a[0][1], a[1][2], a[0][2]};
	}

	// inline Real6x6 stiffnessTensor(Real& shear, Real& bulk){
	// 	return Real6x6{Real6{	2 * shear + bulk, bulk, bulk,	0,	0,	0,},
	// 				   Real6{   bulk, 2 * shear + bulk, bulk,	0,	0,	0,},
	// 				   Real6{	bulk, bulk, 2 * shear + bulk,	0,	0,	0,},
	// 				   Real6{	0,	0,	0, bulk,	0,	0,},
	// 				   Real6{	0,	0,	0,	0, bulk,	0,},
	// 				   Real6{	0,	0,	0,	0,	0, bulk,}};
	// }

	inline void setDims(Real3x3& a,Uint dims){
		if (dims == 2){
										   a[0][2] = 0.0; 
										   a[1][2] = 0.0; 
			a[2][0] = 0.0;  a[2][1] = 0.0; a[2][2] = 1.0;
		} else if (dims == 1){
							a[0][1] = 0.0; a[0][2] = 0.0;
			a[1][0] = 0.0;  a[1][1] = 1.0; a[1][2] = 0.0;
			a[2][0] = 0.0;  a[2][1] = 0.0; a[2][2] = 1.0;
		}
	}

	inline void setDims2(Real3x3& a,Uint dims){
		if (dims == 2){
			a[2][2] = 1.0;
		} else if (dims == 1){
			a[1][1] = 1.0;
			a[2][2] = 1.0;
		}
	}

	inline Real3x3 inv3x3( Real3x3& a, Uint dims){
		Real3x3 b;
		Real a11 = a[0][0], a12 = a[0][1], a13 = a[0][2];
		Real a21 = a[1][0], a22 = a[1][1], a23 = a[1][2];
		Real a31 = a[2][0], a32 = a[2][1], a33 = a[2][2];
		Real d = a11 * (a33 * a22 - a32 * a23) - a12 * (a21 * a33 - a23 * a31) + a13 * (a21 * a32 - a22 * a31);
		b[0][0] =  (a33 * a22 - a32 * a23)/d,  b[0][1] = - (a33 * a12 - a32 * a13)/d, b[0][2] =  (a23 * a12 - a22 * a13)/d;
		b[1][0] = -(a33 * a21 - a31 * a23)/d,  b[1][1] =   (a33 * a11 - a31 * a13)/d, b[1][2] = -(a23 * a11 - a21 * a13)/d;
		b[2][0] =  (a32 * a21 - a31 * a22)/d,  b[2][1] = - (a32 * a11 - a31 * a12)/d, b[2][2] =  (a22 * a11 - a21 * a12)/d;
		return b;
	}

	inline Real det3x3( Real3x3& a){
		return    a[0][0] * (a[2][2] * a[1][1] - a[2][1] * a[1][2]) 
			    - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0]) 
			    + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
	}

	inline Real trace( Real3x3& a){	
		return a[0][0] + a[1][1] + a[2][2];
	}

	inline Real3x3 identity(){
		return Real3x3{Real3{1.0,0.0,0.0},Real3{0.0,1.0,0.0},Real3{0.0,0.0,1.0}};
	}

 	inline Real3x3 tensorProduct(Real3& a, Real3& b){
 		Real3x3 c;

 		c[0][0] = a[0]*b[0]; c[0][1] = a[0]*b[1]; c[0][2] = a[0]*b[2];
 		c[1][0] = a[1]*b[0]; c[1][1] = a[1]*b[1]; c[1][2] = a[1]*b[2];
 		c[2][0] = a[2]*b[0]; c[2][1] = a[2]*b[1]; c[2][2] = a[2]*b[2];


		// if (c[0][0] <= 1.0E-10 && c[0][0] >= -1.0E-10 ) c[0][0] = 1.0;
		// if (c[1][1] <= 1.0E-10 && c[1][1] >= -1.0E-10 ) c[1][1] = 1.0;
		// if (c[2][2] <= 1.0E-10 && c[2][2] >= -1.0E-10 ) c[2][2] = 1.0;

 		return c;
 	}

	inline Real3x3x3 tensorProduct(Real3x3 a, Real3 b){
		Real3x3x3 c;
		for(int i=0;i<3;i++) for(int j=0;j<3;j++) for(int k=0;k<3;k++)
			c[i][j][k] = a[i][j] * b[k];
		return c;
	}

	inline void checkSingularity(Real3x3& L_i){

		if (abs(L_i[0][0]) < 1.0E-10 &&
				abs(L_i[1][0]) < 1.0E-10 && abs(L_i[0][1]) < 1.0E-10 &&
				abs(L_i[2][0]) < 1.0E-10 && abs(L_i[0][2]) < 1.0E-10){
			L_i[0][0] = 1.0; L_i[0][1] = 0; L_i[0][2] = 0;
			L_i[1][0] = 0; 
			L_i[2][0] = 0; 
		}
		if (abs(L_i[1][1]) < 1.0E-10 &&
				abs(L_i[1][2]) < 1.0E-10 && abs(L_i[2][1]) < 1.0E-10 &&
				abs(L_i[0][1]) < 1.0E-10 && abs(L_i[1][0]) < 1.0E-10){
						   L_i[0][1] = 0; 
			L_i[1][0] = 0; L_i[1][1] = 1.0; L_i[1][2] = 0;
						   L_i[2][1] = 0; 
			
		}
		if (abs(L_i[2][2]) < 1.0E-10 &&
				abs(L_i[2][0]) < 1.0E-10 && abs(L_i[2][1]) < 1.0E-10 &&
				abs(L_i[0][2]) < 1.0E-10 && abs(L_i[1][2]) < 1.0E-10){
										  L_i[0][2] = 0;
										  L_i[1][2] = 0;
			L_i[2][0] = 0; L_i[2][1] = 0; L_i[2][2] = 1.0;
			
		}

	}


	inline Real W_Wendland_1D_2h(Real dist, Real smoothingLength){
		Real h = smoothingLength * 0.5;
		Real q = (dist / h);
		if (q > 2.0){
			return 0.0;
		} else{
			return (0.625 / h) * (1.0 - 0.5 * q) *  (1.0 - 0.5 * q) * (1.0 - 0.5 * q) * (1.0 + 1.5 * q);
	   }
	}

 	inline Real W_Wendland_2D_2h(Real dist, Real smoothingLength){
 		Real h = smoothingLength * 0.5;
 		Real q = (dist / h);
 		if (q > 2.0){
 			return 0.0;
 		} else{
 			Real foo = (1.0 - 0.5 * q);
			Real bar = foo * foo;
			Real res = (0.55704230082163367519109317180380 / (h * h)) * bar * bar * (2.0 * q + 1.0);
			return res;
		}
 	}

	inline Real W_CubicSpline_2D_h(Real dist, Real smoothingLength){
 		Real q = (dist / smoothingLength);
		Real bar = (1.0-q);
 		if (q > 1.0){
 			return  0;
 		} else if (q > 0.5){
			return 3.6378272706718933890 * bar * bar * bar;
		} else{
			Real foo = (0.5 - q);
			return 3.6378272706718933890 * (bar * bar * bar - 4.0 * foo * foo * foo);
		}
	}


	inline Real W_Wendland6_2D_h(Real dist, Real smoothingLength){
		Real h = smoothingLength / 2.0;
		Real q = (dist / h);
		if(q > 2.0){
			return 0.0;			
		} else{
			Real foo = (1.0 - q / 2.0);
			Real bar = foo * foo;
			return (0.88672039722627401357 / (h * h)) * (bar * bar * bar * bar) * (1.0 + 4.0 * q + 6.25 * q * q + 4.0 * q * q * q);
		}
	}

	inline Real3 gW_Wendland6_2D_h(Real dist, Real3 dir, Real smoothingLength){

		Real h = smoothingLength / 2.0;
		if(dist / h > 2.0){
			return Real3{0,0,0};			
		}

		Real foo = (2.0 * h - dist);
		Real foo7 = foo * foo * foo * foo * foo * foo * foo;

		Real h13 = h * h * h * h;
			 h13 = h13 * h13 * h13 * h;

		Real3 res = mult(
							- 0.019050633534158230760 * foo7 * dist * (2.0 * h * h + 7.0 * h * dist + 8.0 * dist * dist) / h13,
							  dir
						);

		return res;
	}


 	inline Real W_Wendland_3D_2h(Real dist, Real smoothingLength){
 		Real h = smoothingLength * 0.5;
 		Real q = (dist / h);
 		if (q > 2.0){
 			return 0.0;
 		} else{
 			Real foo = (1.0 - 0.5 * q);
			Real bar = foo * foo;
 			return (0.41778172561622525639 / (h * h * h)) * bar * bar * (2.0 * q + 1.0);
		}
 	}

	
	//  (1-q/2)^4 * (2q + 1)


 	inline Real W_Wendland_2D(Real dist, Real smoothingLength){
 		Real q = dist / smoothingLength;
 		if (q > 1.0){
 			return 0.0;
 		} else{
 			Real foo = (1.0 - q);
			Real bar = foo * foo;
			// return ((4.456338406572776) / (smoothingLength*smoothingLength)) * 0.5 * ((1.0-q)*(1.0-q)*(1.0-q)*(1.0-q)) * (1.0 + 4.0 * q);
	 		return (2.228169203286388) / ( smoothingLength * smoothingLength ) * bar * bar * (1.0 + 4.0 * q);
		}
	 }
	 
	inline Real W_QuinticSpline_3D(Real dist, Real smoothingLength){
		Real r = dist / smoothingLength;
		if ( r > 1.0 ) return 0;
		else if ( 0 < r <= 0.3333333333333333 ) {
			Real a = (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r);
			Real b = (-6.0) * (0.6666666666666666 - r) * (0.6666666666666666 - r) * (0.6666666666666666 - r) * (0.6666666666666666 - r) * (0.6666666666666666 - r);
			Real c = (15.0) * (0.3333333333333333 - r) * (0.3333333333333333 - r) * (0.3333333333333333 - r) * (0.3333333333333333 - r) * (0.3333333333333333 - r);
			return 17.403593027098754966 * (a + b + c);
		} else if ( 0.3333333333333333 < r <= 0.6666666666666666 ) {
			Real a = (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r);
			Real b = (-6.0) * (0.6666666666666666 - r) * (0.6666666666666666 - r) * (0.6666666666666666 - r) * (0.6666666666666666 - r) * (0.6666666666666666 - r);
			return 17.403593027098754966 * (a + b);
		} else{
			Real a = (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r);
			return 17.403593027098754966 * (a);
		}
	}

	inline Real3 gW_QuinticSpline_3D(Real dist, std::array<Real,3> dir, Real smoothingLength){
		Real r = dist / smoothingLength;
		if ( r > 1.0 ) return Real3{0.0,0.0,0.0};
		else if ( 0 < r < 0.3333333333333333) {
			Real a = (-5.0) * (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r);
			Real b = (30.0) * (0.6666666666666666 - r) * (0.6666666666666666 - r) * (0.6666666666666666 - r) * (0.6666666666666666 - r);
			Real c = (-75.0) * (0.3333333333333333 - r) * (0.3333333333333333 - r) * (0.3333333333333333 - r) * (0.3333333333333333 - r);
			return mult((17.403593027098754966 / smoothingLength) * (a + b + c), dir);
		} else if ( 0.3333333333333333 < r <= 0.6666666666666666 ) {
			Real a = (-5.0) * (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r);
			Real b = (30.0) * (0.6666666666666666 - r) * (0.6666666666666666 - r) * (0.6666666666666666 - r) * (0.6666666666666666 - r);			
			return mult((17.403593027098754966 / smoothingLength) * (a + b), dir);
		} else{
			Real a = (-5.0) * (1.0 - r) * (1.0 - r) * (1.0 - r) * (1.0 - r);
			return mult((17.403593027098754966 / smoothingLength) * (a), dir);		
		}
	}

 	inline Real3 gW_Wendland_1D_2h(Real dist, std::array<Real,3> dir, Real smoothingLength){
		Real h = smoothingLength * 0.5;
		Real q = (dist / h);
		if (q > 2.0){
			return Real3{0.0,0.0,0.0};
		} else{
			Real3 res = mult( (-0.625 / (h*h)) * (2.25 + 0.75 * q) * (1.0 - 0.5 * q) * (1.0 - 0.5 * q) , dir);
			return res;
	   }
	}

	inline Real3 gW_Wendland_1D_2h_Nondim(Real dist, std::array<Real,3> dir, Real smoothingLength){
		Real h = smoothingLength * 0.5;
		Real q = (dist / h);
		if (q > 2.0){
			return Real3{0.0,0.0,0.0};
		} else{
			Real3 res = mult( (2.25 + 0.75 * q) * (1.0 - 0.5 * q) * (1.0 - 0.5 * q) , dir);
			return res;
	   }
	}	

	inline Real3 gW_Wendland_3D_2h(Real dist, std::array<Real,3> dir, Real smoothingLength){
 		Real h = smoothingLength * 0.5;
 		Real q = (dist / h);
 		if (q > 2.0){
 			return Real3{0.0,0.0,0.0};
 		} else{
 			Real foo = 1.0 - 0.5 * q;
 			// This is the "2h" version.
 			// (Referred from the DualSphysics Guidbook, \frac{alpha_D}{h} * (-5 * q) * (1-\frac{1}{2}q)^3)
	 		return mult(((0.41778172561622525639)/(h*h*h*h))*(-5.0 * q)*foo*foo*foo, dir);
		}
 	}

	 inline Real3 gW_Wendland_3D_2h_Nondim(Real dist, std::array<Real,3> dir, Real smoothingLength){
		Real h = smoothingLength * 0.5;
		Real q = (dist / h);
		if (q > 2.0){
			return Real3{0.0,0.0,0.0};
		} else{
			Real foo = 1.0 - 0.5 * q;
			// This is the "2h" version.
			// (Referred from the DualSphysics Guidbook, \frac{alpha_D}{h} * (-5 * q) * (1-\frac{1}{2}q)^3)
			return mult((-5.0 * q)*foo*foo*foo, dir);
	   }
	}


 	inline Real3 gW_Wendland_2D_2h(Real dist, std::array<Real,3> dir, Real smoothingLength){
 		Real h = smoothingLength * 0.5;
 		Real q = (dist / h);
 		if (q > 2.0){
 			return Real3{0.0,0.0,0.0};
 		} else{
 			Real foo = 1.0 - 0.5 * q;
 			// This is the "2h" version.
			 // (Referred from the DualSphysics Guidbook, \frac{alpha_D}{h} * (-5 * q) * (1-\frac{1}{2}q)^3)
			 
	 		return mult(((7.0/(4.0*M_PI))/(h*h*h))*(-5.0 * q)*foo*foo*foo, dir);
		}
 	}


 	inline Real3 gW_Wendland_2D_2h_Nondim(Real dist, std::array<Real,3> dir, Real smoothingLength){
		Real h = smoothingLength * 0.5;
		Real q = (dist / h);
		if (q > 2.0){
			return Real3{0.0,0.0,0.0};
		} else{
			Real foo = 1.0 - 0.5 * q;			
			return mult((-5.0 * q)*foo*foo*foo, dir);
	   }
	}	 

 	inline Real3 gW_Wendland_2D(Real dist, std::array<Real,3> dir, Real smoothingLength){
 		Real q = dist / smoothingLength;
 		if (q > 1.0){
 			return Real3{0.0,0.0,0.0};
 		} else{
	 		return mult( ( (4.456338406572776) / ( smoothingLength * smoothingLength * smoothingLength ) ) * (-10.0 * q) * (1.0 - q) * (1.0 - q) * (1.0 - q), dir);
		}
 	}


 	inline Real delta_SPH(Real  rho_i,      Real rho_j,
 						  Real  vol_j,      Real delta,
 						  Real  soundSpeed, Real smoothingLength,
 						  Real  mass_j,     Real  dist,
 						  Real3 relpos,    	Real3 gWij,
						  Real3 densGrad_i, Real3 densGrad_j){
 		// const Real psi = 2.0 * ( - rho_j + rho_i) - dot(add(densGrad_j,densGrad_i),relpos);
 		// const Real qoo = dot(relpos, gWij) * vol_j * soundSpeed * smoothingLength * delta / (dist * dist);
		// return psi * qoo;
		return 2.0 * delta * smoothingLength * soundSpeed * (rho_i - rho_j) * (dot(relpos, gWij) / (dist * dist)) * (mass_j / rho_j);
 	}

 	// inline Real temperatureDependentViscosity(Real mu, Real T){return 0;}
 	inline Real equal (Real T){
 		return T;
 	}

 	inline Real fixedViscosity(Real mu, Real T){return mu;}

	inline Real3x3 tensorViscosity_ij(Real Cl, Real Cq, Real smoothingLength, Real soundSpeed, 
								   Real3x3 velGrad_i, Real3x3 velGrad_j){

		Real3x3 velGrad_ij = sub(velGrad_i, velGrad_j); 
		return transpose(  
						   add( mult(-Cl * smoothingLength * soundSpeed, transpose(velGrad_ij)),
								mult( Cq * smoothingLength * smoothingLength, mult(velGrad_ij, velGrad_ij))
				              )
						);
	} 

	inline Real scalarViscosity_ij(Real Cl, Real Cq, Real smoothingLength, Real soundSpeed, 
								   Real3 relPos, Real3 relVel){
		Real3 eta = mult(1./(smoothingLength), relPos);
		Real mu = (1.0 / (length2(eta) + 0.001 * smoothingLength)) * dot(relVel, eta) ;
		return (-Cl * soundSpeed * mu + Cq * mu * mu);
	}

 	inline Real linearEOS(Real T, Real T0, Real rho, Real rho0, Real c, Real alpha){
 		return c * c * ( (rho - rho0) + alpha * (T - T0) );
 	}

 	inline Real taitEOS(Real T, Real T0, Real rho, Real rho0, Real c, Real alpha){
 		return (1.0/7.0) * rho0 * c * c * (std::pow((rho/rho0),7) - 1.0);
 	}

	inline Real3 iif_acc(Real sij, Real3 reldir, Real dist, Real h, Real m_i, Real m_j){
		return mult( - (sij / m_i) * (m_i + m_j) * cos(1.5 * PI * dist / (3.0 * h)), reldir);
	}

 	inline Real3 viscosity_acc_ij_Shao(Real mj, Real3 relpos, Real3 gWij, Real rho_i, Real rho_j, Real dist, Real3 relvel, Real mu){

 		return mult(( 4.0 * mj * mu * dot(relpos, gWij) / ((rho_i + rho_j) * (dist * dist + EPSL_SMALL )) ), relvel );
 	}


	inline Real3 extraStress_acc_ij(Real mj, Real3 gWij, Real rho_i, Real rho_j, Real3x3 tau_i, Real3x3 tau_j){
		return 
		mult(
			add(mult(1./(rho_i*rho_i),tau_i),mult(1./(rho_j*rho_j),tau_j)),
			mult(mj, gWij)
		);
	}

 	inline Real3 pressure_acc_ij_1(Real mj, Real3 gWij, Real P_i, Real P_j, Real rho_i, Real rho_j, Real vol_j){
 		return mult( - mj * (P_i + P_j) / (rho_i * rho_j), gWij);
 	}

 	inline Real3 pressure_acc_ij_2(Real mj, Real3 gWij, Real P_i, Real P_j, Real rho_i, Real rho_j, Real vol_j){
 		return mult( - (1.0/rho_i) * (P_i + P_j) * vol_j, gWij);
 	}

 	inline Real3 pressure_acc_ij_3(Real mj, Real3 gWij, Real P_i, Real P_j, Real rho_i, Real rho_j, Real vol_j){
 		return mult( - mj * ( P_i / (rho_i * rho_i) + P_j / (rho_j * rho_j) ), gWij);
 	}

 	inline Real constantConductivity(std::string type, Real T){
 		return 1.0;
 	}

	inline Real3 gradient(Real quantityi, Real quantityj, Real3 gWij, Real volj){
		return mult((quantityj-quantityi)*volj, gWij);
	}

 	inline Real cleary(Real Ti, Real Tj, Real mj, Real rho_i, Real rho_j, Real ki, Real kj, Real3 relpos, Real dist, Real3 gWij){
 		return 4.0 * mj * ki * kj * (Ti - Tj) * dot(relpos, gWij) / (rho_i * rho_j * (ki + kj) * (dist * dist));
	 }

	inline Real inconsistentHeatTransfer(
				  					     Real Ti, Real Tj, Real mj, Real rho_i, Real rho_j,
									     Real3 relpos, Real3 reldir,
									     Real dist, Real3 gWij, Real vol_j){
		return 2.0 * ((Ti-Tj) / dist) * dot(reldir,gWij) * vol_j;		
		// return dot(sub(gradT_j,gradT_i),gWij) * vol_j;		
	}
	
	inline Real consistentHeatTransfer(Real3x3 L2_i, Real3 gradT_i,
						   			   Real Ti, Real Tj, Real mj, Real rho_i, Real rho_j, Real ki, Real kj,
						   			   Real3 relpos, Real3 reldir,
						   			   Real dist, Real3 gWij, Real vol_j){
		Real3x3 eij_gWij = tensorProduct(reldir,gWij);

		// return 4.0 * mj * ki * kj * (Ti - Tj) * dot(relpos, gWij) / (rho_i * rho_j * (ki + kj) * (dist * dist));
		// return 4.0 * (ki * kj / (ki + kj)) / (rho_i) * doubleDot(L2_i,eij_gWij) * ((Ti - Tj) / (dist) - (ki * kj / (ki + kj)) * dot(reldir,gradT_i)) * vol_j;
		// return (4.0 / (rho_i)) * doubleDot(L2_i,eij_gWij) * (ki * kj / (ki + kj)) * ((Ti - Tj) / (dist) - dot(reldir,gradT_i)) * vol_j;
		
		return 2.0 * doubleDot(L2_i,eij_gWij) * ((Ti - Tj) / (dist) - dot(reldir,gradT_i)) * vol_j;

		// return 2.0 * ((Ti-Tj) / dist) * dot(reldir,gWij) * vol_j;
	}

 	inline Real3 gravity_acc(){
 		return Real3{0.0,-10.0,0.0};
 	}

 	struct SymTensor3{
 		Real data[3][3][3];
 		SymTensor3(){
 			for(int i=0;i<3;i++)
 				for(int j=0;j<3;j++)
 					for(int k=0;k<3;k++)
 						data[i][j][k] = 0.0;
 		}
 		void add(Uint k, Uint m, Uint n, Real value){
 			data[k][m][n] += value;
 		}
 		void mult(Uint k, Uint m, Uint n, Real value){
 			data[k][m][n] = data[k][m][n] * value;
 		}
 		inline Real operator()(Uint k, Uint m, Uint n){
 			return data[k][m][n];
 		}
 	};
	 
	struct Tensor3{
		Real data[3][3][3];

		Tensor3(){
			for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
			for(int k=0;k<3;k++)
				data[i][j][k] = 0.0;
		}

		Real& operator()(Uint i, Uint j, Uint k){
			return data[i][j][k];
		}
				
	};		

	struct Tensor4{
		Real data[3][3][3][3];

		Tensor4(){
			for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
			for(int k=0;k<3;k++)
			for(int l=0;l<3;l++)
				data[i][j][k][l] = 0.0;
		}

		Real& operator()(Uint i, Uint j, Uint k, Uint l){
			return data[i][j][k][l];
		}

		void add(Tensor4& b){
			for(int i=0;i<3;i++)
			for(int j=0;j<3;j++)
			for(int k=0;k<3;k++)
			for(int l=0;l<3;l++)
				data[i][j][k][l] = data[i][j][k][l] + b(i,j,k,l);			
		}
		
	};		 


// mnk,kl,lop->mnop
	inline void contract_323_to_4(Tensor3& a, Real3x3& b, Tensor3&c, Tensor4& d){
		for(int m=0;m<3;m++) for(int n=0;n<3;n++) for(int k=0;k<3;k++) for(int l=0;l<3;l++) for(int o=0;o<3;o++) for(int p=0;p<3;p++)
			d(m,n,o,p) = d(m,n,o,p) + a(m,n,k) * b[k][l] * c(l,o,p);
	}


}

#endif
