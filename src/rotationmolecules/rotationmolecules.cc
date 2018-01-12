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

/***************************************************************************************/  
#define PI 3.141592653589

int main (int argc, char *argv[]) {
   
	ScreenUtils scrut;
			
	OutputAlignProgram output;
	
	if(argc > 1){


		vector<Atom> molecule_A;
		
		ReadXYZFile reader;

		bool statusAllData_molecule_A = reader.getValuesFromFile(argv[1],molecule_A);

		if(statusAllData_molecule_A){

			VectorAndMatrixOperations matrixOP;

			vector<Atom> molecule_B_operate;

			int decision = 0;
			cout << "What kind of operation want?" << endl;
			cout << "Rotation tap 1, Invertion tap 2" << endl;
			cin >> decision;

			switch(decision){
				case 1:{
					double phi; 
					double theta; 
					double psi; 
					vector<double> angles(3,0.0);

					cout << "Angles for rotation " << endl << "phi = " ;
					cin >> phi;
					cout << "theta = " ;
					cin >> theta;
					cout << "psi = ";
					cin >> psi;
					angles[0] = phi = phi * PI / 180.0;
					angles[1] = theta = theta * PI / 180.0;
					angles[2] = psi = psi * PI / 180.0;
					
					molecule_B_operate = matrixOP.rotateMolecule(angles, molecule_A);

					cout << endl << "Coordenates of inital molecule "<< endl;
					output.displayXYZFile(argv[1],molecule_A);

					cout << endl << "Coordenates of molecule after rotation "<< endl;
					reader.sortingAtoms(molecule_B_operate);
					output.displayXYZFile(argv[1],molecule_B_operate);
					string filaname = argv[1];
					filaname += "_operated"; 
					string comment = "Molecule rotated";
					output.saveXYZFile(filaname,comment,molecule_B_operate);
		 
				}break;
				case 2:{

					molecule_B_operate = matrixOP.inversionOfCoordinates(molecule_A);

					cout << endl << "Coordenates of inital molecule "<< endl;
					output.displayXYZFile(argv[1],molecule_A);

					cout << endl << "Coordenates of molecule after inversion "<< endl;
					reader.sortingAtoms(molecule_B_operate);
					output.displayXYZFile(argv[1],molecule_B_operate);
					string filaname = argv[1];
					filaname += "_operated"; 
					string comment = "Molecule rotated";
					output.saveXYZFile(filaname,comment,molecule_B_operate);

				}break;
				default:
					cout << "Dont needed to use this program " << endl; 
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


