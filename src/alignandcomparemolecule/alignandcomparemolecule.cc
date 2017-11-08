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
		cout << "Comparation molecules." << endl;
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
			
			vector<vector<double>> inertiatensor_molecula_A(3,vector<double>(3,0));
			vector<vector<double>> inertiatensor_molecula_B(3,vector<double>(3,0));

			MolecularOperations molecularOP;
			
			vector<Atom> molecule_A_inCM = molecularOP.moveCM2Origin(molecule_A);
			vector<Atom> molecule_B_inCM = molecularOP.moveCM2Origin(molecule_B);

			inertiatensor_molecula_A = molecularOP.inertiaTensor(molecule_A_inCM);
			inertiatensor_molecula_B = molecularOP.inertiaTensor(molecule_B_inCM);
			
			vector<vector<double>> diagmatrix_molecule_A(3,vector<double>(3,0.0));
			vector<vector<double>> eigvectors_molecule_A(3,vector<double>(3,0.0));
			vector<double> eigvalues_molecule_A(3,0.0);

			vector<vector<double>> diagmatrix_molecule_B(3,vector<double>(3,0.0));
			vector<vector<double>> eigvectors_molecule_B(3,vector<double>(3,0.0));
			vector<double> eigvalues_molecule_B(3,0.0);

			VectorAndMatrixOperations matrixOP;

			matrixOP.eigenVectorValues(inertiatensor_molecula_A,diagmatrix_molecule_A,eigvectors_molecule_A,eigvalues_molecule_A);
			matrixOP.eigenVectorValues(inertiatensor_molecula_B,diagmatrix_molecule_B,eigvectors_molecule_B,eigvalues_molecule_B);

			cout << endl << " Inertia Tensor - Molecule A" << setw(49);
			cout << " Inertia Tensor - Molecule B" << endl;

			for(int i=0;i<3;++i){ 
				cout << " | " << setw(15) << inertiatensor_molecula_A[0][i] << setw(15) << inertiatensor_molecula_A[1][i] << setw(13) << inertiatensor_molecula_A[2][i] << " | ";
				cout << " | " << setw(15) << inertiatensor_molecula_B[0][i] << setw(15) << inertiatensor_molecula_B[1][i] << setw(13) << inertiatensor_molecula_B[2][i] << " | " << endl;
			}
			
			cout << endl << " Diagonalization - Inertia Tensor - Molecule A" << setw(49);
			cout  << " Diagonalization - Inertia Tensor - Molecule B" << endl;

			for(int i=0;i<3;++i){
				cout << " | " << setw(15) << diagmatrix_molecule_A[0][i] << setw(15) << diagmatrix_molecule_A[1][i] << setw(13) << diagmatrix_molecule_A[2][i] << " | ";
				cout << " | " << setw(15) << diagmatrix_molecule_B[0][i] << setw(15) << diagmatrix_molecule_B[1][i] << setw(13) << diagmatrix_molecule_B[2][i] << " | " << endl;
			}
			
			cout << endl << " EingenValues - Inertia Tensor - Molecule A" << setw(49);
			cout << " EingenValues - Inertia Tensor - Molecule B" << endl;

			cout << " | "  << setw(15)<< eigvalues_molecule_A[0]  << setw(15) << eigvalues_molecule_A[1]  << setw(13)<< eigvalues_molecule_A[2] << " | ";
			cout << " | " << eigvalues_molecule_B[0]  << setw(15) << eigvalues_molecule_B[1] << setw(13)<<  eigvalues_molecule_B[2] << " | " << endl << endl;

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


