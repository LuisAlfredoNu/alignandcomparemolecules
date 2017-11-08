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
		vector<double> rotationOperation(double, vector<double>);
		double getAngleBetween2Vectors(vector<double>,vector<double>);
		void eigenVectorValues(vector<vector<double>> initialmatrix,vector<vector<double>> &diagmatrix, vector<vector<double>>& eigvectors,vector<double>& eigvalues);
	/***************************************************************************************/
	/***************************************************************************************/ 
	private:
		//vector<vector<double>> transposeMatrix(vector<vector<double>> matrix);

	/***************************************************************************************/ 
	/***************************************************************************************/ 

	protected:
};
#endif // _VECTORMATRIXOPERATIONS_H_
