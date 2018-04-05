#include <cstdlib>
#include <iostream>
using std::cout;
using std::endl;
using std::left;
using std::cerr;
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

	if(argc > 2){

		output.displayFilesNames(argv[1],argv[2]);

		vector<Atom> molecule_A;
		vector<Atom> molecule_B;

		ReadXYZFile reader;

		bool statusAllData_molecule_A = reader.getValuesFromFile(argv[1],molecule_A);
		bool statusAllData_molecule_B = reader.getValuesFromFile(argv[2],molecule_B);

		if(statusAllData_molecule_A && statusAllData_molecule_B){

			MolecularOperations molecularOP;

			if(molecularOP.haveSameTypeNumAtoms(molecule_A,molecule_B)){

				vector<vector<double>> inertiatensor_molecula_A(3,vector<double>(3,0));
				vector<vector<double>> inertiatensor_molecula_B(3,vector<double>(3,0));

				vector<Atom> molecule_A_inCM = molecularOP.moveCM2Origin(molecule_A);
				vector<Atom> molecule_B_inCM = molecularOP.moveCM2Origin(molecule_B);
				/*
					string filaname_moleculeA = argv[1];
					filaname_moleculeA = filaname_moleculeA.substr(0,filaname_moleculeA.size()-4);
					filaname_moleculeA += "_in_CM.xyz";
					string comment = "Move to Center of mass";

					output.saveXYZFile(filaname_moleculeA,"Molecule A",molecule_A_inCM);

					string filaname_moleculeB = argv[2];
					filaname_moleculeB = filaname_moleculeB.substr(0,filaname_moleculeB.size()-4);
					filaname_moleculeB += "_in_CM.xyz";

					output.saveXYZFile(filaname_moleculeB,"Molecule B",molecule_B_inCM);
					*/
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

				string title = "Inertia Tensor";
				output.displayDualMatrix(title,inertiatensor_molecula_A,inertiatensor_molecula_B);

				title = "EingenVectors - Inertia Tensor";
				output.displayDualMatrix(title,eigvectors_molecule_A,eigvectors_molecule_B);

				title = "EingenValues - Inertia Tensor";
				output.displayDualMatrix(title,eigvalues_molecule_A,eigvalues_molecule_B);

				if(matrixOP.compareEigenValues(eigvalues_molecule_A,eigvalues_molecule_B)){

					output.displayItsTheSame();

					vector<Atom> molecule_A_align = matrixOP.rotateMolecule(eigvectors_molecule_A,molecule_A_inCM);
					vector<Atom> molecule_B_align = matrixOP.rotateMolecule(eigvectors_molecule_B,molecule_B_inCM);

					reader.sortingAtoms(molecule_A_align);
					reader.sortingAtoms(molecule_B_align);

					if(matrixOP.compareCoordinates(molecule_A_align,molecule_B_align)){
						cout << "Equal" << endl;
						//	cout << "The both molecules are the same 00 "<< endl;
					}else{

						vector<double> angles (3,0.0);
						bool is_same_after_rotation_in_Z = false;
						vector<Atom> molecule_B_align_second_rotation;

						while(! is_same_after_rotation_in_Z){

							angles[2] += 90.0;
							molecule_B_align_second_rotation = matrixOP.rotateMolecule(angles,molecule_B_align);

							if(matrixOP.compareCoordinates(molecule_A_align,molecule_B_align_second_rotation) || angles[2] == 270.0)
								is_same_after_rotation_in_Z = true;
						}
						molecule_B_align = molecule_B_align_second_rotation;
						/*
							vector<vector<double>> change2_A_basis = matrixOP.changeBasisEigenVec(eigvectors_molecule_A,eigvectors_molecule_B);
							vector<Atom> molecule_B_align_second_rotation = matrixOP.rotateMolecule2(change2_A_basis,molecule_B_inCM);

							title = "EingenVectors - Change Basis";
							output.displayDualMatrix(title,eigvectors_molecule_A,change2_A_basis);

							inertiatensor_molecula_A = molecularOP.inertiaTensor(molecule_A_align);
							inertiatensor_molecula_B = molecularOP.inertiaTensor(molecule_B_align_second_rotation);

							matrixOP.eigenVectorValues(inertiatensor_molecula_A,diagmatrix_molecule_A,eigvectors_molecule_A,eigvalues_molecule_A);
							matrixOP.eigenVectorValues(inertiatensor_molecula_B,diagmatrix_molecule_B,eigvectors_molecule_B,eigvalues_molecule_B);

							string title = "Inertia Tensor";
							output.displayDualMatrix(title,inertiatensor_molecula_A,inertiatensor_molecula_B);

							title = "EingenVectors - Inertia Tensor";
							output.displayDualMatrix(title,eigvectors_molecule_A,eigvectors_molecule_B);

							molecule_B_align = matrixOP.rotateMolecule(eigvectors_molecule_B,molecule_B_align_second_rotation);
						//molecule_B_align = molecule_B_align_second_rotation;
						*/
						reader.sortingAtoms(molecule_A_align);
						reader.sortingAtoms(molecule_B_align);

						if(matrixOP.compareCoordinates(molecule_A_align,molecule_B_align)){
							cout << "Equal 222" << endl;
							//	cout << "The both molecules are the same 11 "<< endl;
						}else{
							cout << "Stereoisomer" << endl;
							//	cout << "The both molecules are isomers " << endl;
						}
					}

				}else{
					//output.displayDifferentMolecule();
					cout << "Different" << endl;
				}
				cout << endl << "Coordenates of molecule A" << endl;
				//output.saveXYZFile(argv[1],"Molecule A",molecule_A_align);
				output.displayXYZFile(argv[1],molecule_A_inCM);

				cout << endl << "Coordenates of molecule B" << endl;
				//output.saveXYZFile(argv[2],"Molecule B",molecule_B_align);
				output.displayXYZFile(argv[2],molecule_B_inCM);
				
				return EXIT_SUCCESS;
			}else{
				cout << endl;
				scrut.PrintScrStarLine();
				cout << "Different" << endl;
				//scrut.DisplayErrorMessage(" The Molecules are differents ");
				scrut.PrintScrStarLine();
				cout << endl;
				return EXIT_SUCCESS;
			}
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


