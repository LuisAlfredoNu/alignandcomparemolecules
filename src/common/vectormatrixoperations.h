/***************************************************************************************/ 
/*
 * 
 * Class for transform diferent types of arrays 
 * Methods
 * 	getAngleBetween2Vectors -> return the angle between two vector with origin in the same point 
 *
 */
/***************************************************************************************/  
#ifndef _VECTORMATRIXOPERATIONS_H_
#define _VECTORMATRIXOPERATIONS_H_

#include <vector>
using std::vector;
/***************************************************************************************/ 
#include "vectormatrixoperations.h"
/***************************************************************************************/ 
/***************************************************************************************/ 

class VectorAndMatrixOperations{
	/***************************************************************************************/ 
	public:
		VectorAndMatrixOperations();
	/***************************************************************************************/ 
		void eigenVectorValues(vector<vector<double>> initialmatrix,vector<vector<double>> &diagmatrix, vector<vector<double>>& eigvectors,vector<double>& eigvalues);
		bool compareEigenValues(vector<double>, vector<double>);
      bool proveEigenValDegeneracy(vector<double>);

		vector<Atom> rotateMolecule(vector<double> angles,vector<Atom> molecule);
		vector<Atom> rotateMolecule(vector<vector<double>>,vector<Atom> molecule);
		vector<Atom> rotateMolecule2(vector<vector<double>>,vector<Atom> molecule);
		vector<Atom> inversionOfCoordinates(vector<Atom> molecule);
		bool compareCoordinates(vector<Atom>,vector<Atom>);
		bool permutationBequalA(vector<Atom>,vector<Atom> &);
		
		vector<vector<double>> changeBasisEigenVec(vector<vector<double>>,vector<vector<double>>);
		
		vector<double> rotationOperationOverZ(double, vector<double>);
		vector<double> rotationOperationOverY(double, vector<double>);
		vector<double> rotationOperationOverX(double, vector<double>);
		double getAngleBetween2Vectors(vector<double>,vector<double>);
		vector<double> alignEigenVectors4Angles(vector<vector<double>>);

		vector<vector<double>> incrementLengthVector(double,vector<vector<double>>);
		double RMS4Comparations(vector<Atom> molecule_A, vector<Atom> molecule_B);
	/***************************************************************************************/
	/***************************************************************************************/ 
	private:
		double dotProduct(vector<double>,vector<double>);
		vector<double> crossProduct(vector<double>,vector<double>);
		vector<double> anglesEuler(int,vector<vector<double>>);
		vector<vector<double>> rotationEuler(vector<double>, vector<vector<double>> ); 

	/***************************************************************************************/ 
	/***************************************************************************************/ 

	protected:
};
#endif // _VECTORMATRIXOPERATIONS_H_
