#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::left;
using std::cerr;
#include <iomanip>
using std::setw;
#include <vector>
using std::vector;
/***************************************************************************************/ 
#include "screenutils.h"
#include "readxyzfile.h"
#include "atomsinmolecule.h"
#include "molecularoperations.h"
#include "vectormatrixoperations.h"

int main (int argc, char *argv[]) {
   
	ScreenUtils scrut;
	
	if(argc > 2){

		cout << endl;
		scrut.PrintScrStarLine();
		scrut.SetScrGreenBoldFont();
		cout << "Rotation molecules." << endl;
		scrut.SetScrNormalFont();
		scrut.PrintScrStarLine();
		cout << endl << " Molecule A: "<< argv[1] << setw(35);
		cout << " Molecule B: "<< argv[2] << endl;

		vector<Atom> molecule_A;
		vector<Atom> molecule_B;
		
		ReadXYZFile reader;

		bool statusAllData_molecule_A = reader.getValuesFromFile(argv[1],molecule_A);
		bool statusAllData_molecule_B = reader.getValuesFromFile(argv[2],molecule_B);

		if(statusAllData_molecule_A && statusAllData_molecule_B){
			
			
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
		cout << "Example: rotationmolecules file_A.xyz file_B.xyz  "<< endl;
		scrut.PrintScrStarLine();
		cout << endl;
		return EXIT_FAILURE;
	}
}


