#include <cstdlib>
#include <iostream>
using std::cin;
using std::cout;
using std::endl;
using std::left;
using std::cerr;
#include <iomanip>
using std::setw;
#include <vector>
using std::vector;
#include <string>
using std::string;
/***************************************************************************************/ 
#include "screenutils.h"
#include "readxyzfile.h"
#include "atomsinmolecule.h"
#include "molecularoperations.h"
#include "vectormatrixoperations.h"
#include "output.h"

int main (int argc, char *argv[]) {
   
	ScreenUtils scrut;
			
	OutputAlignProgram output;
	
	if(argc > 1){


		vector<Atom> molecule_A;
		
		ReadXYZFile reader;

		bool statusAllData_molecule_A = reader.getValuesFromFile(argv[1],molecule_A);

		if(statusAllData_molecule_A){

			VectorAndMatrixOperations matrixOP;

			double phi; 
			double theta; 
			double psi; 
			
			cout << "Angles for rotation " << endl << "phi = " ;
			cin >> phi;
			cout << "theta = " ;
			cin >> theta;
			cout << "psi = ";
			cin >> psi;
			phi = phi * 3.14159 / 180.0;
			theta = theta * 3.14159 / 180.0;
			psi = psi * 3.14159 / 180.0;
			vector<double> angles = {phi,theta,psi};

			vector<Atom> molecule_B_rotate = matrixOP.rotateMolecule(angles, molecule_A);

			cout << endl << "Coordenates of inital molecule "<< endl;
			output.displayXYZFile(argv[1],molecule_A);
			
			cout << endl << "Coordenates of molecule after rotation "<< endl;
			reader.sortingAtoms(molecule_B_rotate);
			output.displayXYZFile(argv[1],molecule_B_rotate);
			string filaname = argv[1];
			filaname += "_rotated"; 
			string comment = "Molecule rotated";
			output.saveXYZFile(filaname,comment,molecule_B_rotate);
			
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


