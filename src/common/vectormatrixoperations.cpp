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
#include <math.h>
#define PI 	3.141592653589793
/***************************************************************************************/ 
#include"eig2-4.h"
#include"vectormatrixoperations.h"
/***************************************************************************************/  
/***************************************************************************************/  

VectorAndMatrixOperations::VectorAndMatrixOperations(){ }
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
double VectorAndMatrixOperations::getAngleBetween2Vectors(vector<double> vector01, vector<double> vector02){
	double dotproduct = 0.0;
	for(int xyz=0;xyz<3;xyz++) dotproduct += vector01[xyz] * vector02[xyz];
	
	double magnitudVec01 = 0.0, magnitudVec02 = 0.0;
	for(int xyz=0;xyz<3;xyz++){
		magnitudVec01 += vector01[xyz] * vector01[xyz]; 
		magnitudVec02 += vector02[xyz] * vector02[xyz];
	}
	magnitudVec01 = sqrt(magnitudVec01);
	magnitudVec02 = sqrt(magnitudVec02);

	double angle = 0.0;
	angle = acos(dotproduct / (magnitudVec01*magnitudVec02)) * 180.0 / PI;

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
	vector<double> magnitudA (3,0.0);
	vector<double> magnitudB (3,0.0);

	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			magnitudA[i] += eigenVector_moleculeA[i][j] * eigenVector_moleculeA[i][j];
			magnitudB[i] += eigenVector_moleculeB[i][j] * eigenVector_moleculeB[i][j];
		}
		magnitudA[i] = sqrt(magnitudA[i]);
		magnitudB[i] = sqrt(magnitudB[i]);
	}

	int biggestVectorA = -1;
	if(magnitudA[0] > magnitudA[1]){
		if(magnitudA[0] > magnitudA[2]){
			biggestVectorA = 0;
		}else{
			biggestVectorA = 2;
		}
	}else{
		if(magnitudA[1]>magnitudA[2]){
			biggestVectorA = 1;
		}else{
			biggestVectorA = 2;
		}
	}
	int biggestVectorB = -1;
	if(magnitudB[0] > magnitudB[1]){
		if(magnitudB[0] > magnitudB[2]){
			biggestVectorB = 0;
		}else{
			biggestVectorB = 2;
		}
	}else{
		if(magnitudB[1]>magnitudB[2]){
			biggestVectorB = 1;
		}else{
			biggestVectorB = 2;
		}
	}

	vector<double> unitvectorZ (3,0.0);
	unitvectorZ[2] = 1.0;

	double angle_eigVecA = getAngleBetween2Vectors(eigenVector_moleculeA[biggestVectorA], unitvectorZ);
	double angle_eigVecB = getAngleBetween2Vectors(eigenVector_moleculeB[biggestVectorB], unitvectorZ);

	for(int i=0;i<3;i++){
		eigenVector_moleculeA[i] = rotationOperationOverY(angle_eigVecA,eigenVector_moleculeA[i]);
		eigenVector_moleculeB[i] = rotationOperationOverY(angle_eigVecB,eigenVector_moleculeB[i]);
	}

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

