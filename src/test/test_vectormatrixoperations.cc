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


	return EXIT_SUCCESS;
}


