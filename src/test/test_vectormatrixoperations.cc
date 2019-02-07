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
	cout << "********************************************************" << endl;

	vector<double> vector01={5.0, 5.0, 0.0};
	vector<double> vector02={0.0,10.0, 0.0};

	VectorAndMatrixOperations matrixOP;

	double angle = matrixOP.getAngleBetween2Vectors(vector01,vector02);

	cout << endl << "Vector A: " << vector01[0] << " " << vector01[1] << " " << vector01[2] << endl;
	cout << endl << "Vector B: " << vector02[0] << " " << vector02[1] << " " << vector02[2] << endl;

	cout << endl << "Angle between this vector is = " << angle << endl;

	vector<double> rotatedvector (3,0.0);

	rotatedvector = matrixOP.rotationOperationOverZ(45.0, vector01);
	cout << endl << "Vector rotated 45° in the plane x-y in counterwise clock = "<<rotatedvector[0]<<" , "<<rotatedvector[1]<<" , "<<rotatedvector[2] << endl;
	
	rotatedvector = matrixOP.rotationOperationOverY(45.0, vector01);
	cout << endl << "Vector rotated 45° in the plane x-z in counterwise clock = "<<rotatedvector[0]<<" , "<<rotatedvector[1]<<" , "<<rotatedvector[2] << endl;
	
	rotatedvector = matrixOP.rotationOperationOverX(45.0, vector01);
	cout << endl << "Vector rotated 45° in the plane y-z in counterwise clock = "<<rotatedvector[0]<<" , "<<rotatedvector[1]<<" , "<<rotatedvector[2] << endl;

	vector<vector<double>> initialmatrix (3,vector<double>(3,0.0));
	vector<vector<double>> diagmatrix(3,vector<double>(3,0.0));
	vector<vector<double>> eigvectors(3,vector<double>(3,0.0));
	vector<vector<double>> transpose_eigvectors(3,vector<double>(3,0.0));
	vector<double> eigvalues(3,0.0);

  initialmatrix[0][0]=446379	;	initialmatrix[0][1]=-38006.1;	initialmatrix[0][2]=94002.7; 
  initialmatrix[1][0]=-38006.1;	initialmatrix[1][1]=405182	;	initialmatrix[1][2]=155591 ;
  initialmatrix[2][0]=94002.7	;	initialmatrix[2][1]=155591	;	initialmatrix[2][2]=87697.2; 


  matrixOP.eigenVectorValues(initialmatrix,diagmatrix,eigvectors,eigvalues);
  
  cout << endl << " Inertia Tensor - Matrix" << endl;
  for(int i=0;i<3;++i) cout << " | " << initialmatrix[0][i] << "\t--\t" << initialmatrix[1][i]<< "\t--\t" << initialmatrix[2][i] << " | " << endl;

  cout << endl << " EingenVectores - Inertia Tensor - Matrix" << endl;
  for(int i=0;i<3;++i) cout << " | " << eigvectors[0][i] << "\t--\t" << eigvectors[1][i]<< "\t--\t" << eigvectors[2][i] << " | " << endl    ;
  
  cout << endl << " EingenValues - Inertia Tensor - Matrix" << endl;
  cout << " | " << eigvalues[0] << "\t--\t" << eigvalues[1]<< "\t--\t" << eigvalues[2] << " | " << endl;

  cout << endl << " Diagonalization - Inertia Tensor - Matrix" << endl;
  for(int i=0;i<3;++i) cout << " | " << diagmatrix[0][i] << "\t--\t" << diagmatrix[1][i]<< "\t--\t" << diagmatrix[2][i] << " | " << endl    ;

  bool status = matrixOP.compareEigenValues(eigvalues,eigvalues);

  cout << endl << "Is the same the last eigenVector: " << (status ? "Yes" : "No") << endl << endl; 

	return EXIT_SUCCESS;
}


