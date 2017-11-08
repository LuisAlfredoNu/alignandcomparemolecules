#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <string>
using std::string;
#include <vector>
using std::vector;
#include "vectormatrixoperations.h"

int main (int argc, char *argv[])
{
	cout << endl << "********************************************************" << endl;
	cout << " Testing for VectorAndMatrixOperations Class " << endl;
	cout << "********************************************************" << endl << endl;

	vector<double> vector01={5.0, 5.0, 0.0};
	vector<double> vector02={0.0,10.0, 0.0};

	VectorAndMatrixOperations matrixOP;

	double angle = matrixOP.getAngleBetween2Vectors(vector01,vector02);

	cout << endl << "Angle between this vector is = " << angle << endl;

	vector<double> rotatedvector (3,0.0);
	rotatedvector = matrixOP.rotationOperation(45.0, vector01);

	cout << endl << "Vector rotated 45° in the plane x-y in counterwise clock = "<<rotatedvector[0]<<" , "<<rotatedvector[1]<<" , "<<rotatedvector[2] << endl;

	vector<vector<double>> initialmatrix (3,vector<double>(3,0.0));
	vector<vector<double>> diagmatrix(3,vector<double>(3,0.0));
	vector<vector<double>> eigvectors(3,vector<double>(3,0.0));
	vector<double> eigvalues(3,0.0);

  initialmatrix[0][0]=446379	;	initialmatrix[0][1]=-38006.1;	initialmatrix[0][2]=94002.7; 
  initialmatrix[1][0]=-38006.1;	initialmatrix[1][1]=405182	;	initialmatrix[1][2]=155591 ;
  initialmatrix[2][0]=94002.7	;	initialmatrix[2][1]=155591	;	initialmatrix[2][2]=87697.2; 


  matrixOP.eigenVectorValues(initialmatrix,diagmatrix,eigvectors,eigvalues);

  cout << endl << "Inertia Tensor Diagonalizaded - Matrix" << endl;
  for(int i=0;i<3;++i) cout << " | " << diagmatrix[i][0] << "\t--\t" << diagmatrix[i][1]<< "\t--\t" << diagmatrix[i][2] << " | " << endl    ;



	return EXIT_SUCCESS;
}


