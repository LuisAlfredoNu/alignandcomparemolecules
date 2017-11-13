/***************************************************************************************/ 
/*
 * 
 * Class for transform diferent types of arrays 
 * Methods
 * 	getAngleBetween2Vectors -> return the angle between two vector with origin in the same point 
 *
 */
/***************************************************************************************/  
#include "atomsinmolecule.h"
#ifndef _VECTORMATRIXOPERATIONS_H_
#define _VECTORMATRIXOPERATIONS_H_

#include <vector>
using std::vector;
/***************************************************************************************/ 
#include "atomsinmolecule.h"
/***************************************************************************************/ 
/***************************************************************************************/ 

class VectorAndMatrixOperations{
	/***************************************************************************************/ 
	public:
		VectorAndMatrixOperations();
	/***************************************************************************************/ 
		void eigenVectorValues(vector<vector<double>> initialmatrix,vector<vector<double>> &diagmatrix, vector<vector<double>>& eigvectors,vector<double>& eigvalues);
		bool compareEigenValues(vector<double>, vector<double>);
		void alignEigenVectors(vector<vector<double>> eigenVector_moleculeA,vector<vector<double>> eigenVector_moleculeB);
		vector<double> anglesEuler(int,vector<vector<double>>);
		vector<vector<double>> rotationEuler(vector<double>, vector<vector<double>> ); 
	/***************************************************************************************/
	/***************************************************************************************/ 
	private:
		double dotProduct(vector<double>,vector<double>);
		vector<double> crossProduct(vector<double>,vector<double>);
		vector<double> rotationOperationOverZ(double, vector<double>);
		vector<double> rotationOperationOverY(double, vector<double>);
		vector<double> rotationOperationOverX(double, vector<double>);
		double getAngleBetween2Vectors(vector<double>,vector<double>);
		//vector<vector<double>> transposeMatrix(vector<vector<double>> matrix);

	/***************************************************************************************/ 
	/***************************************************************************************/ 

	protected:
};
#endif // _VECTORMATRIXOPERATIONS_H_
