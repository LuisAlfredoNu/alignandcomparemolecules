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
	cout << "********************************************************" << endl;
	cout << "********************************************************" << endl;
	
	vector<Atom> molecule (5,Atom());

	double x1=1.0, y1=1.0, z1=0.0;
	double x2=-1.0, y2=1.0, z2=0.0;
	double x3=0.0, y3=-1.0, z3=1.0;
	double x4=2.0, y4=-3.0, z4=-0.0;
	double x5=-12.4, y5=-22.9, z5=55.5;
	string nameelemnt01 ("Au");
	int numelement02=4;
	string nameelemnt03 ("H");
	int numelement04=43;
	int numelement05=83;

	molecule[0].setCoordinates(x1,y1,z1);
	molecule[0].setAtomSymbol(nameelemnt01);
	molecule[1].setCoordinates(x2,y2,z2);
	molecule[1].setAtomNumber(numelement02);
	molecule[2].setCoordinates(x3,y3,z3);
	molecule[2].setAtomSymbol(nameelemnt03);
	molecule[3].setCoordinates(x4,y4,z4);
	molecule[3].setAtomNumber(numelement04);
	molecule[4].setCoordinates(x5,y5,z5);
	molecule[4].setAtomNumber(numelement05);

/***************************************************************************************/ 
/***************************************************************************************/ 

	MolecularOperations molecularop;
	vector<double> rmasscenter(3,0);
	vector<vector<double>> inertiatensor(3,vector<double>(3,0));
	vector<Atom> molecule_inCM;

	cout << endl << "********************************************************" << endl;
	cout << " Compute the center of mass" << endl;
	cout << "********************************************************" << endl << endl;

	rmasscenter = molecularop.massCenter(molecule);

	for(int i=0;i<3;++i) cout << "Coordinate of Center of mass = " << rmasscenter[i] << endl;
	
	cout << endl << "********************************************************" << endl;
	cout << " Compute the Inertia Tensor " << endl;
	cout << "********************************************************" << endl << endl;
	
	inertiatensor = molecularop.inertiaTensor(molecule);
	
	cout << endl << "Inertia Tensor - Matrix" << endl;
	for(int i=0;i<3;++i) cout << " | " << inertiatensor[i][0] << "\t--\t" << inertiatensor[i][1]<< "\t--\t" << inertiatensor[i][2] << " | " << endl;

	cout << endl << "********************************************************" << endl;
	cout << " Move the center of mass of the molecule to the origin " << endl;
	cout << "********************************************************" << endl << endl;
	
	molecule_inCM = molecularop.moveCM2Origin(molecule);
	string statusanswer_before, statusanswer_after;

	for(unsigned int i=0;i<molecule_inCM.size();i++){
		cout << "Atom = " << i << endl;
		for (int ii=0;ii<3;ii++) cout << "Coordinate before = " << molecule[i].atomCoordinates[ii] << "\t Coordinate after = "<< molecule_inCM[i].atomCoordinates[ii] << endl;
		cout << "Before move ***** Element: "<< molecule[i].atomSymbol <<"   Atomic number="<< molecule[i].atomNumber << "   Atomic weight= "<< molecule[i].atomWeight    << endl;
		cout << "After move  ***** Element: "<< molecule_inCM[i].atomSymbol <<"   Atomic number="<< molecule_inCM[i].atomNumber << "   Atomic weight= "<< molecule_inCM[i].atomWeight    << endl;
		molecule[i].statusData ? statusanswer_before="yes" : statusanswer_before="no";
		molecule_inCM[i].statusData ? statusanswer_after="yes" : statusanswer_after="no";
		cout << " Before All data is fine? " << statusanswer_before << endl << " After All data is fine? " << statusanswer_after << endl;
		cout << "********************************************************" << endl << endl;
	}
	
	rmasscenter = molecularop.massCenter(molecule_inCM);

	cout << " After Center of Mass" << endl;
	for(int i=0;i<3;++i) cout << "Coordinate of Center of mass = " << rmasscenter[i] << endl;

	inertiatensor = molecularop.inertiaTensor(molecule_inCM);
	
	cout << endl<< " After Inertia Tensor" << endl;
	cout << "Inertia Tensor - Matrix" << endl;
	for(int i=0;i<3;++i) cout << " | " << inertiatensor[i][0] << "\t--\t" << inertiatensor[i][1]<< "\t--\t" << inertiatensor[i][2] << " | " << endl;


	return EXIT_SUCCESS;
}


