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

vector<double> VectorAndMatrixOperations::rotationOperation(double angle,vector<double> vector2rotate){
	vector<vector<double>> matrixrotation (3,vector<double> (3,0.0));
	matrixrotation[0][0] = matrixrotation[1][1] = cos(angle * PI / 180.0);
	matrixrotation[0][1] = -sin(angle * PI /180.0);
	matrixrotation[1][0] = sin(angle * PI /180.0);
	matrixrotation[2][0] = matrixrotation[2][1] = matrixrotation[0][2] = matrixrotation[1][2] = 0.0;
	matrixrotation[2][2] = 1.0;

	vector<double> vector_prime(3,0.0);

	for(int i=0;i < 3; i++){
		for(int j=0;j<3;j++){
			vector_prime[i] += matrixrotation[i][j]*vector2rotate[j]; 
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

