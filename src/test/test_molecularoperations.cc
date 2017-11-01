#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <string>
using std::string;
#include <vector>
using std::vector;
#include "atomsinmolecule.h"
#include "molecularoperations.h"

int main (int argc, char *argv[])
{
	cout << endl << "********************************************************" << endl;
	cout << " Testing for MolecularOperations Class " << endl;
	cout << "********************************************************" << endl << endl;
	
	vector<Atom> molecule (4,Atom());

	double x1=1.0, y1=1.0, z1=0.0;
	double x2=-1.0, y2=1.0, z2=0.0;
	double x3=0.0, y3=-1.0, z3=1.0;
	double x4=0.0, y4=-1.0, z4=-1.0;
	string nameelemnt01 ("H");
	int numelement02=1;
	string nameelemnt03 ("H");
	int numelement04=1;

	molecule[0].setCoordinates(x1,y1,z1);
	molecule[0].setAtomSymbol(nameelemnt01);
	molecule[1].setCoordinates(x2,y2,z2);
	molecule[1].setAtomNumber(numelement02);
	molecule[2].setCoordinates(x3,y3,z3);
	molecule[2].setAtomSymbol(nameelemnt03);
	molecule[3].setCoordinates(x4,y4,z4);
	molecule[3].setAtomNumber(numelement04);

/***************************************************************************************/ 
/***************************************************************************************/ 

	MolecularOperations molecularop;
	vector<double> rmasscenter(3,0);
	vector<vector<double>> inertiatensor(3,vector<double>(3,0));

	rmasscenter = molecularop.massCenter(molecule);

	for(int i=0;i<3;++i) cout << "Coordinate of Center of mass = " << rmasscenter[i] << endl;
	
	inertiatensor= molecularop.inertiaTensor(molecule);
	
	cout << endl << "Inertia Tensor - Matrix" << endl;
	for(int i=0;i<3;++i) cout << " | " << inertiatensor[i][0] << "\t--\t" << inertiatensor[i][1]<< "\t--\t" << inertiatensor[i][2] << " | " << endl;


	return EXIT_SUCCESS;
}


