/***************************************************************************************/  
/* Class of OutputAlignProgram       */
/***************************************************************************************/  
#ifndef _OUTPUT_CPP_
#define _OUTOUT_CPP_

#include <iostream>
using std::cout;
using std::endl;
#include <iomanip>
using std::setw;
#include <vector>
using std::vector;
#include <string>
using std::string;
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
	scrut.SetScrGreenBoldFont();
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
/***************************************************************************************/ 
/***************************************************************************************/ 
/***************************************************************************************/ 
#endif // _OUTPUT_CPP_

