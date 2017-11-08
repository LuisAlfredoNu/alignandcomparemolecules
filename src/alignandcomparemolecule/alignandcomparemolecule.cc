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
#include "vectormatrixoperations.h"

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

			vector<vector<double>> initialmatrix_molecule_A (3,vector<double>(3,0.0));
			vector<vector<double>> diagmatrix_molecule_A(3,vector<double>(3,0.0));
			vector<vector<double>> eigvectors_molecule_A(3,vector<double>(3,0.0));
			vector<double> eigvalues_molecule_A(3,0.0);

			vector<vector<double>> initialmatrix_molecule_B (3,vector<double>(3,0.0));
			vector<vector<double>> diagmatrix_molecule_B(3,vector<double>(3,0.0));
			vector<vector<double>> eigvectors_molecule_B(3,vector<double>(3,0.0));
			vector<double> eigvalues_molecule_B(3,0.0);

			VectorAndMatrixOperations matrixOP;

			matrixOP.eigenVectorValues(initialmatrix_molecule_A,diagmatrix_molecule_A,eigvectors_molecule_A,eigvalues_molecule_A);
			matrixOP.eigenVectorValues(initialmatrix_molecule_B,diagmatrix_molecule_B,eigvectors_molecule_B,eigvalues_molecule_B);

			cout << endl << " Diagonalization - Inertia Tensor - Molecule A" << endl;
			for(int i=0;i<3;++i) cout << " | " << diagmatrix_molecule_A[0][i] << "\t--\t" << diagmatrix_molecule_A[1][i]<< "\t--\t" << diagmatrix_molecule_A[2][i] << " | " << endl;
			cout << endl << " EingenValues - Inertia Tensor - Molecule A" << endl;
			cout << " | " << eigvalues_molecule_A[0] << "\t--\t" << eigvalues_molecule_A[1]<< "\t--\t" << eigvalues_molecule_A[2] << " | " << endl;

			cout << endl << " Diagonalization - Inertia Tensor - Molecule B" << endl;
			for(int i=0;i<3;++i) cout << " | " << diagmatrix_molecule_B[0][i] << "\t--\t" << diagmatrix_molecule_B[1][i]<< "\t--\t" << diagmatrix_molecule_B[2][i] << " | " << endl;
			cout << endl << " EingenValues - Inertia Tensor - Molecule B" << endl;
			cout << " | " << eigvalues_molecule_B[0] << "\t--\t" << eigvalues_molecule_B[1]<< "\t--\t" << eigvalues_molecule_B[2] << " | " << endl;


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


