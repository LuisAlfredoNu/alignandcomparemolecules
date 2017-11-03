/***************************************************************************************/  
/* Class of Vector and Matrix Operations       */
/***************************************************************************************/  
#ifndef _VECTORMATRIXOPERATIONS_CPP_
#define _VECTORMATRIXOPERATIONS_CPP_

#include <vector>
using std::vector;
#include <math.h>
#define PI 	3.141592653589793
/***************************************************************************************/ 
#include<vectormatrixoperations.h>
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
/***************************************************************************************/ 
#endif // _VECTORMATRIXOPERATIONS_CPP_

