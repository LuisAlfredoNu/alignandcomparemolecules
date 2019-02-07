#include <vector>
using std::vector;
/***************************************************************************************/ 
#include "screenutils.h"
#include "optflag.h"
#include "readxyzfile.h"
#include "atomsinmolecule.h"
#include "molecularoperations.h"
#include "vectormatrixoperations.h"
#include "output.h"

int main (int argc, char *argv[]) {

	OptFlags optflags;
	optflags.getOptions(argc,argv);

	string filename_molecule_A = argv[argc-2];
	string filename_molecule_B = argv[argc-1];

	OutputAlignProgram output;

	vector<Atom> molecule_A;
	vector<Atom> molecule_B;

   string result;
   double RMSD = 0.0;
   bool degeneracy = false;

	ReadXYZFile reader;
	MolecularOperations molecularOP;
	VectorAndMatrixOperations matrixOP;

	bool statusAllData_molecule_A = reader.getValuesFromFile(filename_molecule_A,molecule_A);
   bool statusAllData_molecule_B = reader.getValuesFromFile(filename_molecule_B,molecule_B);

   output.correctInputData(statusAllData_molecule_A,statusAllData_molecule_B);
	output.displayFilesNames(optflags.quiet_version,filename_molecule_A,filename_molecule_B);

   if(molecularOP.haveSameTypeNumAtoms(molecule_A,molecule_B)){

      vector<vector<double>> inertiatensor_molecula_A(3,vector<double>(3,0));
      vector<vector<double>> inertiatensor_molecula_B(3,vector<double>(3,0));

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

      matrixOP.eigenVectorValues(inertiatensor_molecula_A,diagmatrix_molecule_A,eigvectors_molecule_A,eigvalues_molecule_A);
      matrixOP.eigenVectorValues(inertiatensor_molecula_B,diagmatrix_molecule_B,eigvectors_molecule_B,eigvalues_molecule_B);

      if(optflags.display_inertia_tensor)
         output.displayInertiaTensorEigenVecEigenVal(inertiatensor_molecula_A,inertiatensor_molecula_B,eigvectors_molecule_A,eigvectors_molecule_B,eigvalues_molecule_A,eigvalues_molecule_B);

      if(optflags.display_large_eigenvector)
         output.displayLargeEigenVectors(optflags.vector_increase_length,eigvectors_molecule_A,eigvectors_molecule_B);

      vector<Atom> molecule_A_align;	
      vector<Atom> molecule_B_align;

      if(matrixOP.compareEigenValues(eigvalues_molecule_A,eigvalues_molecule_B)){
         degeneracy = matrixOP.proveEigenValDegeneracy(eigvalues_molecule_A);
         if(!degeneracy){

            molecule_A_align = matrixOP.rotateMolecule(eigvectors_molecule_A,molecule_A_inCM);
            molecule_B_align = matrixOP.rotateMolecule(eigvectors_molecule_B,molecule_B_inCM);

            reader.sortingAtoms(molecule_A_align);
            reader.sortingAtoms(molecule_B_align);

            if(matrixOP.compareCoordinates(molecule_A_align,molecule_B_align)){
               result = "Equal 00 ";
            }else{
               bool find_equal = matrixOP.permutationBequalA(molecule_A_align,molecule_B_align);
               if(find_equal){
                  result = "Equal 11 ";
               }else{
                  result = "Enantiomers ";
               }
            }
            if (optflags.run_big_test){
               inertiatensor_molecula_A = molecularOP.inertiaTensor(molecule_A_align);
               inertiatensor_molecula_B = molecularOP.inertiaTensor(molecule_B_align);
               matrixOP.eigenVectorValues(inertiatensor_molecula_A,diagmatrix_molecule_A,eigvectors_molecule_A,eigvalues_molecule_A);
               matrixOP.eigenVectorValues(inertiatensor_molecula_B,diagmatrix_molecule_B,eigvectors_molecule_B,eigvalues_molecule_B);

               if(optflags.display_inertia_tensor)
                  output.displayInertiaTensorEigenVecEigenVal(inertiatensor_molecula_A,inertiatensor_molecula_B,eigvectors_molecule_A,eigvectors_molecule_B,eigvalues_molecule_A,eigvalues_molecule_B);

               if(optflags.display_large_eigenvector)
                  output.displayLargeEigenVectors(optflags.vector_increase_length,eigvectors_molecule_A,eigvectors_molecule_B);
            }
            if(optflags.display_rms)
               RMSD = matrixOP.RMS4Comparations(molecule_A_align,molecule_B_align);

            if(optflags.display_output_coordenates){
               output.display_booth_XYZFile(filename_molecule_A,filename_molecule_B,molecule_A_align,molecule_B_align);

            }
            if(optflags.save_output_coordenates){
               output.saveXYZFile(filename_molecule_A,"Molecule A",molecule_A_align);
               output.saveXYZFile(filename_molecule_B,"Molecule B",molecule_B_align);
            }
         }else{
            result = "Equal "; 
            if(optflags.display_rms)
               RMSD = matrixOP.RMS4Comparations(molecule_A_inCM,molecule_B_inCM);
         }
      }else{
         result = "Different ";
         if(optflags.display_rms)
            RMSD = matrixOP.RMS4Comparations(molecule_A_inCM,molecule_B_inCM);
      }
   }else{
      result = "Different ";
      if(optflags.display_rms)
         RMSD = matrixOP.RMS4Comparations(molecule_A,molecule_B);
   }
   if(optflags.display_rms)
      output.displayRMSD(RMSD);
   if(degeneracy)
      output.displayDegeneracy();
   output.displayResult(optflags.quiet_version,result);
   
	return EXIT_SUCCESS;
}


