/***************************************************************************************/  
/* Class of OutputAlignProgram       */
/***************************************************************************************/  
#ifndef _OUTPUT_CPP_
#define _OUTOUT_CPP_

#include <iostream>
using std::cout;
using std::endl;
#include <iomanip>
using std::left;
using std::setw;
using std::setprecision;
#include <vector>
using std::vector;
#include <string>
using std::string;
#include <fstream>
using std::ofstream;
/***************************************************************************************/ 
#include "screenutils.h"
#include "output.h"
/***************************************************************************************/  
/***************************************************************************************/  

OutputAlignProgram::OutputAlignProgram(){
}

/***************************************************************************************/ 
void OutputAlignProgram::displayDualMatrix(string title,vector<vector<double>> matrixA,vector<vector<double>> matrixB){

	cout << endl <<" " << title <<" - Molecule A" << setw(title.size()+(-title.size() - 14 + 50));
	cout << title << " - Molecule B" << endl;
	for(unsigned int i=0;i<matrixA[0].size();++i){
		cout << " | " << setw(15) << matrixA[0][i] << setw(15) << matrixA[1][i] << setw(13) << matrixA[2][i] << " | ";
		cout << " | " << setw(15) << matrixB[0][i] << setw(15) << matrixB[1][i] << setw(13) << matrixB[2][i] << " | " << endl;
	}
}
/***************************************************************************************/ 
void OutputAlignProgram::displayDualMatrix(string title,vector<double> matrixA,vector<double> matrixB){

	cout << endl <<" " << title <<" - Molecule A" << setw(title.size()+(-title.size() - 14 + 50));
	cout << title << " - Molecule B" << endl;
		cout << " | " << setw(15) << matrixA[0] << setw(15) << matrixA[1] << setw(13) << matrixA[2] << " | ";
		cout << " | " << setw(15) << matrixB[0] << setw(15) << matrixB[1] << setw(13) << matrixB[2] << " | " << endl;
}
/***************************************************************************************/ 
void OutputAlignProgram::displayFilesNames(string file1,string file2){

	cout << endl;
	scrut.PrintScrStarLine();
	scrut.SetScrYellowBoldFont();
	cout << "Comparation molecules." << endl;
	scrut.SetScrNormalFont();
	scrut.PrintScrStarLine();
	cout << endl << " Molecule A: "<< file1 << setw(file1.size()+(-file1.size() - 11 + 50));
	cout << " Molecule B: "<< file2 << endl;

}
/***************************************************************************************/ 
void OutputAlignProgram::displayItsTheSame(){
	
	scrut.PrintScrStarLine();
	scrut.SetScrGreenBoldFont();
	cout << "Have the same EigenValues" << endl;
	scrut.SetScrNormalFont();
	scrut.PrintScrStarLine();
}
/***************************************************************************************/ 
bool OutputAlignProgram::saveXYZFile(string filename,string commnet, vector<Atom> molecule){

	string ofilename = filename.substr(0,(filename.size()-4));
	ofilename += "_align.xyz";
	ofstream ofil(ofilename.c_str());
	
	if ( !ofil.good() ) {
		cout << "Error while opening the file " << ofilename << endl;
      cout << "OWBCOutputData values not written!" << endl;
		cout << __FILE__ << ", line: " << __LINE__ << endl;
		ofil.close();
		return false;
	}
	ofil << molecule.size() << endl;
	ofil << filename << " ---> " << ofilename << commnet << endl;
	for(unsigned int i=0;i<molecule.size();i++){

		ofil << setw(5) << left << molecule[i].atomSymbol;
		ofil << setw(15)<< left << molecule[i].atomCoordinates[0];
		ofil << setw(15)<< left << molecule[i].atomCoordinates[1];
		ofil << setw(15)<< left << molecule[i].atomCoordinates[2];
		ofil << endl;
	}
	cout << "Output file XYZ " << ofilename << endl;
	return true;
}
/***************************************************************************************/ 
void OutputAlignProgram::displayXYZFile(string filename, vector<Atom> molecule){
	
	cout << endl << "Molecule from file: " << filename << endl;
	cout << molecule.size() << endl;
	cout << endl;
	for(unsigned int i=0;i<molecule.size();i++){

		cout << setw(5) << left << molecule[i].atomSymbol;
		cout << setw(15) << setprecision(5) << molecule[i].atomCoordinates[0];
		cout << setw(15) << setprecision(5) << molecule[i].atomCoordinates[1];
		cout << setw(15) << setprecision(5) << molecule[i].atomCoordinates[2];
		cout << endl;
	}
}
/***************************************************************************************/ 
void OutputAlignProgram::displaySameMolecule(){
	
	scrut.PrintScrStarLine();
	scrut.SetScrGreenBoldFont();
	cout << "Are the same molecule" << endl;
	scrut.SetScrNormalFont();
	scrut.PrintScrStarLine();
}
/***************************************************************************************/ 
void OutputAlignProgram::displayDifferentMolecule(){
	
	scrut.PrintScrStarLine();
	scrut.SetScrRedBoldFont();
	cout << "Are different molecule" << endl;
	scrut.SetScrNormalFont();
	scrut.PrintScrStarLine();
}
/***************************************************************************************/ 
void OutputAlignProgram::displayIsomerMolecule(){
	
	scrut.PrintScrStarLine();
	scrut.SetScrBlueBoldFont();
	cout << "Are Isomers" << endl;
	scrut.SetScrNormalFont();
	scrut.PrintScrStarLine();
}
/***************************************************************************************/ 
/***************************************************************************************/ 
#endif // _OUTPUT_CPP_

