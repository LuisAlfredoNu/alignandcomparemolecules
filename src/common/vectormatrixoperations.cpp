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
	rotation[0][0] = cos(theta * PI / 180.0); 
	rotation[0][1] = -sin(theta * PI /180.0);
	rotation[0][2] = 0.0; 
	rotation[1][0] = sin(theta * PI /180.0);
	rotation[1][1] = cos(theta * PI / 180.0);
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
	rotation[0][0] = cos(phi * PI / 180.0); 
	rotation[0][1] = 0.0;
	rotation[0][2] = sin(phi * PI /180.0);
	rotation[1][0] = 0.0;
	rotation[1][1] = 1.0;
	rotation[1][2] = 0.0;
	rotation[2][0] = -sin(phi * PI /180.0);  
	rotation[2][1] = 0.0; 
	rotation[2][2] = cos(phi * PI / 180.0);
	
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
	rotation[1][1] = cos(psi * PI / 180.0);
	rotation[1][2] = -sin(psi * PI /180.0);
	rotation[2][0] = 0.0; 
	rotation[2][1] = sin(psi * PI /180.0);
	rotation[2][2] = cos(psi * PI / 180.0); 

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
	angle = acos(dotProduct(vectorA,vectorB) / (magnitudVec01*magnitudVec02)) * 180.0 / PI;

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
		transpose_eigvectors[i][j] = array_eigvectors[j][i];
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
void VectorAndMatrixOperations::alignEigenVectors(vector<vector<double>>& eigenVector_moleculeA,vector<vector<double>>& eigenVector_moleculeB){

	vector<double> crossproductA = crossProduct(eigenVector_moleculeA[0],eigenVector_moleculeA[1]);
	vector<double> crossproductB = crossProduct(eigenVector_moleculeB[0],eigenVector_moleculeB[1]);

	cout << "Cross product A v1xv2 = "<< crossproductA[0] << " , " << crossproductA[1] << " , " <<crossproductA[2] << endl;
	cout << "Cross product B v1xv2 = "<< crossproductB[0] << " , " << crossproductB[1] << " , " <<crossproductB[2] << endl;
	
	vector<vector<double>> alignvector (3,vector<double> (3,0.0));
	vector<vector<double>> unitvectorXYZ (3,vector<double> (3,0.0));
	unitvectorXYZ[0][0] = 1.0;
	unitvectorXYZ[1][1] = 1.0;
	unitvectorXYZ[2][2] = 1.0;

	// Align in axis Z
	// Get the projection over x-y axis
	vector<double> component_xy (3,0.0);
	for(int xyz=0;xyz<3;xyz++) {
		component_xy[0] = dotProduct(eigenVector_moleculeA[2], unitvectorXYZ[0]);
		component_xy[1] = dotProduct(eigenVector_moleculeA[2], unitvectorXYZ[1]);
	}

	double angle_eigVecA_xy = getAngleBetween2Vectors(component_xy,unitvectorXYZ[1]);

	alignvector[2] = rotationOperationOverZ(angle_eigVecA_xy,eigenVector_moleculeA[2]);
	
	cout << "Angle of xy A = " << angle_eigVecA_xy << endl;
	cout << endl << " Align Vectors - Molecule A" << endl;
	for(int i=0;i<3;++i){
		cout << " | " << setw(15) << alignvector[0][i] << setw(15) << alignvector[1][i] << setw(13) << alignvector[2][i] << " | " << endl;
	}


	// Align in axis Y
	vector<double> component_yz (3,0.0);
	for(int xyz=0;xyz<3;xyz++) {
		component_yz[1] += eigenVector_moleculeA[2][xyz] * unitvectorXYZ[1][xyz];
		component_yz[2] += eigenVector_moleculeA[2][xyz] * unitvectorXYZ[2][xyz];
	}
	
	double angle_eigVecA_yz = getAngleBetween2Vectors(eigenVector_moleculeA[2],unitvectorXYZ[2]);
	eigenVector_moleculeA[2] = rotationOperationOverX(angle_eigVecA_yz,eigenVector_moleculeA[0]);

}
/*
vector<vector<double>> VectorAndMatrixOperations::transposeMatrix(vector<vector<double>> matrix){
	
	vector<vector<double>> transposematrix (3,vector<double>(3,0.0));

	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){


		}
	}
	return transposematrix;
}
*/

/***************************************************************************************/ 
/***************************************************************************************/ 
#endif // _VECTORMATRIXOPERATIONS_CPP_

