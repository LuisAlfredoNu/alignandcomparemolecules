/***************************************************************************************/  
/* Class of Vector and Matrix Operations       */
/***************************************************************************************/  
#ifndef _VECTORMATRIXOPERATIONS_CPP_
#define _VECTORMATRIXOPERATIONS_CPP_

#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <vector>
using std::vector;
#include <iomanip>
using std::setw;
#include <math.h>
#define PI 	3.141592653589793
/***************************************************************************************/ 
#include"eig2-4.h"
#include"atomsinmolecule.h"
#include"vectormatrixoperations.h"
/***************************************************************************************/  
/***************************************************************************************/  

VectorAndMatrixOperations::VectorAndMatrixOperations(){ }
/***************************************************************************************/  
double VectorAndMatrixOperations::dotProduct(vector<double> vectorA, vector<double> vectorB){
	double resultdotproduct = 0.0;
	for(int xyz=0;xyz<3;xyz++) resultdotproduct += vectorA[xyz] * vectorB[xyz];
	return resultdotproduct;
}
vector<double> VectorAndMatrixOperations::crossProduct(vector<double> vectorA, vector<double> vectorB){
	vector<double> crossproduct (3,0.0);
	crossproduct[0] = vectorA[1]*vectorB[2] - vectorA[2]*vectorB[1]; 
	crossproduct[1] = vectorA[2]*vectorB[0] - vectorA[0]*vectorB[2]; 
	crossproduct[2] = vectorA[0]*vectorB[1] - vectorA[1]*vectorB[0]; 
	return crossproduct;
}
/*
 * For rotation operation 
 * All rotation is counter clockwise
 */
vector<double> VectorAndMatrixOperations::rotationOperationOverZ(double theta,vector<double> vector2rotate){
	vector<vector<double>> rotation (3,vector<double> (3,0.0));
	rotation[0][0] = cos(theta ); 
	rotation[0][1] = sin(theta);
	rotation[0][2] = 0.0; 
	rotation[1][0] = -sin(theta );
	rotation[1][1] = cos(theta );
	rotation[1][2] = 0.0;
	rotation[2][0] = 0.0; 
	rotation[2][1] = 0.0; 
	rotation[2][2] = 1.0;
	
	vector<double> vector_prime(3,0.0);

	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			vector_prime[i] += rotation[i][j]*vector2rotate[j]; 
		}
	}
	return vector_prime;
}
/***************************************************************************************/ 
vector<double> VectorAndMatrixOperations::rotationOperationOverY(double phi,vector<double> vector2rotate){
	vector<vector<double>> rotation (3,vector<double> (3,0.0));
	rotation[0][0] = cos(phi); 
	rotation[0][1] = 0.0;
	rotation[0][2] = sin(phi);
	rotation[1][0] = 0.0;
	rotation[1][1] = 1.0;
	rotation[1][2] = 0.0;
	rotation[2][0] = -sin(phi);  
	rotation[2][1] = 0.0; 
	rotation[2][2] = cos(phi);
	
	vector<double> vector_prime(3,0.0);

	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			vector_prime[i] += rotation[i][j]*vector2rotate[j]; 
		}
	}
	return vector_prime;
}
/***************************************************************************************/ 
vector<double> VectorAndMatrixOperations::rotationOperationOverX(double psi,vector<double> vector2rotate){
	vector<vector<double>> rotation (3,vector<double> (3,0.0));
	rotation[0][0] = 1.0;
	rotation[0][1] = 0.0;
	rotation[0][2] = 0.0; 
	rotation[1][0] = 0.0; 
	rotation[1][1] = cos(psi );
	rotation[1][2] = sin(psi);
	rotation[2][0] = 0.0; 
	rotation[2][1] = -sin(psi );
	rotation[2][2] = cos(psi ); 

	vector<double> vector_prime(3,0.0);

	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			vector_prime[i] += rotation[i][j]*vector2rotate[j]; 
		}
	}
	return vector_prime;
}
/***************************************************************************************/ 
double VectorAndMatrixOperations::getAngleBetween2Vectors(vector<double> vectorA, vector<double> vectorB){
	double magnitudVec01 = 0.0, magnitudVec02 = 0.0;
	for(int xyz=0;xyz<3;xyz++){
		magnitudVec01 += vectorA[xyz] * vectorA[xyz]; 
		magnitudVec02 += vectorB[xyz] * vectorB[xyz];
	}
	magnitudVec01 = sqrt(magnitudVec01);
	magnitudVec02 = sqrt(magnitudVec02);

	double angle = 0.0;
	angle = acos(dotProduct(vectorA,vectorB) / (magnitudVec01*magnitudVec02)) ;

	return angle;
}
/***************************************************************************************/ 
void VectorAndMatrixOperations::eigenVectorValues(vector<vector<double>> initialmatrix,vector<vector<double>> &diagmatrix, vector<vector<double>>& eigvectors,vector<double>& eigvalues){ 
	double array_initialmatrix[3][3];
	double array_eigvectors[3][3];
	double array_eigvalues[3]; 

	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			array_initialmatrix[i][j] = initialmatrix[i][j];
		}
	}

	eigen_decomposition3(array_initialmatrix,array_eigvectors,array_eigvalues);

	vector<vector<double>> transpose_eigvectors(3,vector<double>(3,0.0));

	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
		eigvectors[i][j] = array_eigvectors[i][j];
		transpose_eigvectors[i][j] = array_eigvectors[i][j];
		}
		eigvalues[i] = array_eigvalues[i];
	}

	vector<vector<double>> pre_diagmatrix(3,vector<double>(3,0.0));
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			for(int k=0;k<3;k++){
				pre_diagmatrix[i][j] += transpose_eigvectors[i][k] * initialmatrix[k][j];  
			}
		}
	}
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			for(int k=0;k<3;k++){
				diagmatrix[i][j] += pre_diagmatrix[i][k]* eigvectors[k][j];  
			}
		}
	}
}
/***************************************************************************************/ 
bool VectorAndMatrixOperations::compareEigenValues(vector<double> eigenValues_moleculeA, vector<double> eigenValues_moleculeB){
	bool is_equal = true;
	double epsilon = 0.5;
	double diffvalues = 0.0;

	for(int i=0;i<3;i++){
		diffvalues = eigenValues_moleculeA[i] - eigenValues_moleculeB[i];
		diffvalues = abs(diffvalues);
		if(diffvalues > epsilon) is_equal = false;
	}
	return is_equal;
}
/***************************************************************************************/ 
vector<double> VectorAndMatrixOperations::alignEigenVectors4Angles(vector<vector<double>> eigenVector_molecule){
	
	vector<vector<double>> unitvectorXYZ (3,vector<double> (3,0.0));
	unitvectorXYZ[0][0] = 1.0;
	unitvectorXYZ[1][1] = 1.0;
	unitvectorXYZ[2][2] = 1.0;
	
	double phi = 0.0;
	double theta = 0.0;
	double psi = 0.0;

	vector<vector<double>> alignvector (3,vector<double>(3,0.0)); 
	vector<vector<double>> tmp_alignvector (3,vector<double>(3,0.0)); 
	// Align the axis Z

	phi = atan2(eigenVector_molecule[2][1],eigenVector_molecule[2][2]);
	vector<double> projection_xy_A = {0.0, eigenVector_molecule[2][1], eigenVector_molecule[2][2] }; 
	phi = getAngleBetween2Vectors(projection_xy_A,unitvectorXYZ[2]);
	for(int i=0;i<3;i++)
		tmp_alignvector[i] = rotationOperationOverX(phi,eigenVector_molecule[i]);
	if(abs(unitvectorXYZ[2][1]-tmp_alignvector[2][1])>1e-8){
		phi = -phi;
		for(int i=0;i<3;i++)
			tmp_alignvector[i] = rotationOperationOverX(phi,eigenVector_molecule[i]);
	}
	alignvector = tmp_alignvector;
	cout << " Phi = " << (phi * 180.0 / PI);

	theta = getAngleBetween2Vectors(alignvector[2],unitvectorXYZ[2]);
	for(int i=0;i<3;i++)
		tmp_alignvector[i] = rotationOperationOverY(theta,alignvector[i]);
		
	if(abs(dotProduct(tmp_alignvector[2],unitvectorXYZ[2]) - 1.0)>1e-8){
		theta = -theta;
		for(int i=0;i<3;i++)
			tmp_alignvector[i] = rotationOperationOverY(theta,alignvector[i]);
	}
	cout << " Theta = " << (theta * 180.0 / PI);
	alignvector = tmp_alignvector;

	psi = getAngleBetween2Vectors(alignvector[1],unitvectorXYZ[1]);
	for(int i=0;i<3;i++)
		tmp_alignvector[i] = rotationOperationOverZ(psi,alignvector[i]);
		
	if(abs(dotProduct(tmp_alignvector[1],unitvectorXYZ[1]) - 1.0)>1e-8){
		psi = -psi;
		for(int i=0;i<3;i++)
			tmp_alignvector[i] = rotationOperationOverZ(psi,alignvector[i]);
	}
	cout << " Psi = " << (psi * 180.0 / PI) << endl;
	alignvector = tmp_alignvector;

	cout << endl << " Align Vectors - Molecule " << endl;
	for(int i=0;i<3;++i){
		cout << " | " << setw(15) << alignvector[0][i] << setw(15) << alignvector[1][i] << setw(13) << alignvector[2][i] << " | " << endl;
	}

	vector<double> angles (3,0.0);
	angles[0] = phi;
	angles[1] = theta;
	angles[2] = psi;

	return angles;
/*
	vector<vector<double>> alignvectorB (3,vector<double>(3,0.0)); 
	vector<vector<double>> tmp_alignvectorB (3,vector<double>(3,0.0)); 
	// Align the axis Z

	phi = atan2(eigenVector_moleculeB[2][1],eigenVector_moleculeB[2][2]);
	vector<double> projection_xy = {0.0, eigenVector_moleculeB[2][1], eigenVector_moleculeB[2][2] }; 
	phi = getAngleBetween2Vectors(projection_xy,unitvectorXYZ[2]);
	for(int i=0;i<3;i++)
		tmp_alignvectorB[i] = rotationOperationOverX(phi,eigenVector_moleculeB[i]);
	if(abs(unitvectorXYZ[2][1]-tmp_alignvectorB[2][1])>1e-8){
		phi = -phi;
		for(int i=0;i<3;i++)
			tmp_alignvectorB[i] = rotationOperationOverX(phi,eigenVector_moleculeB[i]);
	}
	alignvectorB = tmp_alignvectorB;
	cout << " Phi = " << (phi * 180.0 / PI);

	theta = getAngleBetween2Vectors(alignvectorB[2],unitvectorXYZ[2]);
	for(int i=0;i<3;i++)
		tmp_alignvectorB[i] = rotationOperationOverY(theta,alignvectorB[i]);
		
	if(abs(dotProduct(tmp_alignvectorB[2],unitvectorXYZ[2]) - 1.0)>1e-8){	
		theta = -theta;
		for(int i=0;i<3;i++)
			tmp_alignvectorB[i] = rotationOperationOverY(theta,alignvectorB[i]);
	}
	cout << " Theta = " << (theta * 180.0 / PI);
	alignvectorB = tmp_alignvectorB;

	psi = getAngleBetween2Vectors(alignvectorB[1],alignvectorA[1]);
	for(int i=0;i<3;i++)
		tmp_alignvectorB[i] = rotationOperationOverZ(psi,alignvectorB[i]);
		
	if(abs(dotProduct(tmp_alignvectorB[1],unitvectorXYZ[1]) - 1.0)>1e-8){
		psi = -psi;
		for(int i=0;i<3;i++)
			tmp_alignvectorB[i] = rotationOperationOverZ(psi,alignvectorB[i]);
	}
	cout << " Psi = " << (psi * 180.0 / PI) << endl;
	alignvectorB = tmp_alignvectorB;
	
	cout << endl << " Align Vectors - Molecule B" << endl;
	for(int i=0;i<3;++i){
		cout << " | " << setw(15) << alignvectorB[0][i] << setw(15) << alignvectorB[1][i] << setw(13) << alignvectorB[2][i] << " | " << endl;
	}
*/
}
/***************************************************************************************/ 
vector<Atom> VectorAndMatrixOperations::rotateMolecule(vector<double> angles,vector<Atom> molecule){
	vector<Atom> molecule_aligned = molecule;
	for(unsigned int i=0;i < molecule.size();i++){

		molecule_aligned[i].setCoordinates(rotationOperationOverX(angles[0],molecule[i].atomCoordinates));
		molecule_aligned[i].setCoordinates(rotationOperationOverY(angles[1],molecule_aligned[i].atomCoordinates));
		molecule_aligned[i].setCoordinates(rotationOperationOverZ(angles[2],molecule_aligned[i].atomCoordinates));
	}

	return molecule_aligned;
}
/***************************************************************************************/  
vector<Atom> VectorAndMatrixOperations::rotateMolecule(vector<vector<double>> matrixrotation,vector<Atom> molecule){
	vector<Atom> molecule_aligned = molecule;
	vector<vector<double>> matrixrotation_tras (3,vector<double> (3,0.0));
	for(int i=0;i<3;i++){
		matrixrotation_tras[0][i] = matrixrotation[i][0];
		matrixrotation_tras[1][i] = matrixrotation[i][1];
		matrixrotation_tras[2][i] = matrixrotation[i][2];
	}
	for(unsigned int i=0;i < molecule.size();i++){
		vector<double> coordinates (3,0.0);
		vector<double> coordinatesB (3,0.0);
		for(int j=0; j<3;j++){
			coordinates[0] += matrixrotation_tras[0][j] * molecule[i].atomCoordinates[j];
			coordinates[1] += matrixrotation_tras[1][j] * molecule[i].atomCoordinates[j];
			coordinates[2] += matrixrotation_tras[2][j] * molecule[i].atomCoordinates[j];
		}
		/*
		for(int j=0; j<3;j++){
			coordinatesB[0] += coordinates[j] * matrixrotation[j][0] ;
			coordinatesB[1] += coordinates[j] * matrixrotation[j][1] ;
			coordinatesB[2] += coordinates[j] * matrixrotation[j][2] ;
		}
		*/
		molecule_aligned[i].setCoordinates(coordinates); 
	}
	return molecule_aligned;
}
/***************************************************************************************/ 
vector<vector<double>> VectorAndMatrixOperations::changeBasisEigenVec(vector<vector<double>> basisA,vector<vector<double>> basisB){
	
	vector<vector<double>> changed_basis(3,vector<double>(3,0.0));

	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			for(int k=0;k<3;k++){
				changed_basis[i][j] += basisA[i][k] * basisB[k][j]; 
			}
		}
	}

	return changed_basis;
}
/***************************************************************************************/ 
vector<double> VectorAndMatrixOperations::anglesEuler(int numangle,vector<vector<double>> eigenVector_moleculeA){
	vector<double> angles (3,0.0);
	vector<vector<double>> unitvectorXYZ (3,vector<double> (3,0.0));
	unitvectorXYZ[0][0] = 1.0;
	unitvectorXYZ[1][1] = 1.0;
	unitvectorXYZ[2][2] = 1.0;
	
	double phi = 0.0;
	double theta = 0.0;
	double psi = 0.0;

	if(numangle == 0){
		phi = atan2(eigenVector_moleculeA[0][1],eigenVector_moleculeA[0][0]) - atan2(unitvectorXYZ[0][1],unitvectorXYZ[0][0]) ;
	}
	if(numangle == 1){
		theta = atan2(eigenVector_moleculeA[2][2],eigenVector_moleculeA[2][0]) - atan2(unitvectorXYZ[2][2], unitvectorXYZ[2][0]) ;
	}
	if(numangle == 2){
	}
	if(numangle == 3){
		phi = atan2(eigenVector_moleculeA[0][1],eigenVector_moleculeA[0][0]) - atan2(unitvectorXYZ[0][1],unitvectorXYZ[0][0]) ;
		theta = atan2(eigenVector_moleculeA[2][2],eigenVector_moleculeA[2][0]) - atan2(unitvectorXYZ[2][2], unitvectorXYZ[2][0]) ;
		psi = atan2(eigenVector_moleculeA[0][2],eigenVector_moleculeA[0][0]) - atan2(unitvectorXYZ[2][1],unitvectorXYZ[2][0]) ;
	}

	cout << " Phi = " << (phi * 360 /(2.0*PI));
	cout << " Tetha = " << (theta * 180 /PI);
	cout << " Psi = " << (psi * 180 /PI) << endl;
	
	angles[0] = phi;
	angles[1] = theta;
	angles[2] = psi;
	/*	
	double tmp_angle1 = atan2(component_xy[0][0],component_xy[0][1]) - atan2(unitvectorXYZ[0][0], unitvectorXYZ[0][0]);
	angles[0] = 101.3004 * (1.0*PI) / 180.0;
	angles[1] = -24.25 * (1.0*PI) / 180.0;
	angles[2] = 0.0 * (1.0*PI) / 180.0;
	tmp_angle1 *= 180.0 / (1.0 * PI);

	if(tmp_angle1 < 0.0) tmp_angle1 += 360.0;

	double tmp_angle2 = atan2(component_xy[1][0],component_xy[1][1])  - atan2(unitvectorXYZ[1][0], unitvectorXYZ[1][1]);
	tmp_angle2 *= 360.0 / (2.0 * PI);

	if(tmp_angle2 < 0.0) tmp_angle2 += 360.0;
	
	cout << "phi atan2 = " << tmp_angle1 <<" - " << tmp_angle2 << endl;
	cout << "phi = " << angleVec1_2_compXY <<" - " << angleVec2_2_compXY << endl;
	*/
	return angles;
}
/***************************************************************************************/ 
vector<vector<double>> VectorAndMatrixOperations::rotationEuler(vector<double> angleseuler,vector<vector<double>> eigenVector_moleculeA){
	vector<vector<double >> matrixrotation (3,vector<double>(3,0.0));

	matrixrotation[0][0] =  cos(angleseuler[2]) * cos(angleseuler[0]) - cos(angleseuler[1]) * sin(angleseuler[0]) * sin(angleseuler[2]);
	matrixrotation[1][0] = -sin(angleseuler[2]) * cos(angleseuler[0]) - cos(angleseuler[1]) * sin(angleseuler[0]) * cos(angleseuler[2]);
	matrixrotation[2][0] =  sin(angleseuler[1]) * sin(angleseuler[0]);

	matrixrotation[0][1] =  cos(angleseuler[2]) * sin(angleseuler[0]) + cos(angleseuler[1]) * cos(angleseuler[0]) * sin(angleseuler[2]);
	matrixrotation[1][1] = -sin(angleseuler[2]) * sin(angleseuler[0]) + cos(angleseuler[1]) * cos(angleseuler[0]) * cos(angleseuler[2]);
	matrixrotation[2][1] = -sin(angleseuler[1]) * cos(angleseuler[0]);

	matrixrotation[0][2] =  sin(angleseuler[2]) * sin(angleseuler[1]);
	matrixrotation[1][2] =  cos(angleseuler[2]) * sin(angleseuler[1]);
	matrixrotation[2][2] =  cos(angleseuler[1]);


	vector<vector<double>> afterrotation_vec (3,vector<double>(3,0.0));
	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			for(int k=0;k<3;k++){
				afterrotation_vec[i][j] += matrixrotation[j][k]*eigenVector_moleculeA[i][k]; 
			}
		}
	}


	return afterrotation_vec;
}
/***************************************************************************************/ 
vector<Atom> VectorAndMatrixOperations::inversionOfCoordinates(vector<Atom> molecule){

	vector<vector<double>> matrixinvertion (3,vector<double>(3,0.0));
	matrixinvertion[0][0] = -1;
	matrixinvertion[1][1] = -1;
	matrixinvertion[2][2] = -1;

	vector<Atom> molecule_inverse = molecule;
	
	for(unsigned int i=0;i < molecule.size();i++){
		vector<double> coordinates (3,0.0);
		for(int j=0; j<3;j++){
			coordinates[0] += matrixinvertion[0][j] * molecule[i].atomCoordinates[j];
			coordinates[1] += matrixinvertion[1][j] * molecule[i].atomCoordinates[j];
			coordinates[2] += matrixinvertion[2][j] * molecule[i].atomCoordinates[j];
		}
		molecule_inverse[i].setCoordinates(coordinates); 
	}
	return molecule_inverse;
}
/***************************************************************************************/ 
bool VectorAndMatrixOperations::compareCoordinates(vector<Atom> molecule_A, vector<Atom> molecule_B){

	bool is_equal = true;
	unsigned int maxsize = molecule_A.size();
	unsigned int i = 0;

	while(is_equal && i < maxsize ){

		if(abs(molecule_A[i].atomCoordinates[0] - molecule_B[i].atomCoordinates[0]) > 0.001) is_equal = false;
		if(abs(molecule_A[i].atomCoordinates[1] - molecule_B[i].atomCoordinates[1]) > 0.001) is_equal = false;
		if(abs(molecule_A[i].atomCoordinates[2] - molecule_B[i].atomCoordinates[2]) > 0.001) is_equal = false;

		i++;
	}

	return is_equal;
}
/***************************************************************************************/ 
/***************************************************************************************/ 
#endif // _VECTORMATRIXOPERATIONS_CPP_

