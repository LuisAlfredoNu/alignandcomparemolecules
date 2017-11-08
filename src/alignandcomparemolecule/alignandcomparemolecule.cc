#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::cerr;
#include <vector>
using std::vector;
/***************************************************************************************/ 
#include "screenutils.h"
#include "readxyzfile.h"
#include "atomsinmolecule.h"
#include "molecularoperations.h"

int main (int argc, char *argv[]) {
   
	ScreenUtils scrut;
	
	if(argc > 2){
		cout << "File for molecule A:  "<<argv[1]<<endl;
		cout << "File for molecule B:  "<<argv[2]<<endl;

		vector<Atom> molecule_A;
		vector<Atom> molecule_B;
		
		ReadXYZFile reader;

		bool statusAllData_molecule_A = reader.getValuesFromFile(argv[1],molecule_A);
		bool statusAllData_molecule_B = reader.getValuesFromFile(argv[2],molecule_B);

		if(statusAllData_molecule_A && statusAllData_molecule_B){
			
			vector<vector<double>> inertiatensor_molecula_A(3,vector<double>(3,0));
			vector<vector<double>> inertiatensor_molecula_B(3,vector<double>(3,0));

			MolecularOperations molecularOP;
			
			vector<Atom> molecule_A_inCM = molecularOP.moveCM2Origin(molecule_A);
			vector<Atom> molecule_B_inCM = molecularOP.moveCM2Origin(molecule_B);

			inertiatensor_molecula_A = molecularOP.inertiaTensor(molecule_A_inCM);
			inertiatensor_molecula_B = molecularOP.inertiaTensor(molecule_B_inCM);



			return EXIT_SUCCESS;
		}else{
			cout << endl;
			scrut.PrintScrStarLine();
			scrut.DisplayErrorMessage(" Problems to read input files ");
			scrut.PrintScrStarLine();
			cout << endl;
			return EXIT_FAILURE;
		}
	}else{
		cout << endl;
		scrut.PrintScrStarLine();
		scrut.DisplayErrorMessage("Dont input files in arguments ");
		cout << "Example: alignandcomparemolecule file_A.xyz file_B.xyz  "<< endl;
		scrut.PrintScrStarLine();
		cout << endl;
		return EXIT_FAILURE;
	}
}


