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

			double phi = 45.0;
			double theta = 30.0;
			double psi = 12.0;
			phi = phi * 3.14159 / 180.0;
			theta = theta * 3.14159 / 180.0;
			psi = psi * 3.14159 / 180.0;
			vector<double> angles = {phi,theta,psi};

			vector<Atom> molecule_B_rotate = matrixOP.rotateMolecule(angles, molecule_A);

			cout << endl << "Coordenates of inital molecule "<< endl;
			cout << molecule_A.size() << endl << endl;
			for(unsigned int i=0;i<molecule_A.size();i++){
				cout << molecule_A[i].atomSymbol << setw(13);
				cout << molecule_A[i].atomCoordinates[0] <<setw(15) ;
				cout << molecule_A[i].atomCoordinates[1] <<setw(15) ;
				cout << molecule_A[i].atomCoordinates[2] ;
				cout << endl;
			}
			cout << endl << "Coordenates of molecule after rotation "<< endl;
			cout << "\tAngles rotation over Z axis = " << phi << endl;
			cout << "\tAngles rotation over Y axis = " << theta << endl;
			cout << "\tAngles rotation over X axis = " << psi << endl;
			cout << molecule_B_rotate.size() << endl;
			cout << argv[1] << endl;
			for(unsigned int i=0;i<molecule_B_rotate.size();i++){
				cout << molecule_B_rotate[i].atomSymbol << setw(15);
				cout << molecule_B_rotate[i].atomCoordinates[0] << setw(15);
				cout << molecule_B_rotate[i].atomCoordinates[1] << setw(15);
				cout << molecule_B_rotate[i].atomCoordinates[2];
				cout << endl;
			}
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


