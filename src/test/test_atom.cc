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

int main (int argc, char *argv[])
{
	vector<Atom> atoms (3,Atom());

	double x1=1.5, y1=2.0, z1=0.5;
	double x2=3.0, y2=1.0, z2=1.5;
	double x3=1.2, y3=3.4, z3=5.6;
	string nameelemnt01 ("H");
	string nameelemnt02 ("C");
	string nameelemnt03 ("Ne");

	atoms[0].setCoordinates(x1,y1,z1);
	atoms[0].setTypeElement(nameelemnt01);
	atoms[1].setCoordinates(x2,y2,z2);
	atoms[1].setTypeElement(nameelemnt02);
	atoms[2].setCoordinates(x3,y3,z3);
	atoms[2].setTypeElement(nameelemnt03);
	for(int i=0;i<3;i++)
		cout << "Coordinates of atom "<< i <<" = " << atoms[i].getXCoordinate() <<","<< atoms[i].getYCoordinate() <<","<< atoms[i].getZCoordinate() <<"  Element: "<< atoms[i].getTypeElement() << endl;


	return EXIT_SUCCESS;
}


