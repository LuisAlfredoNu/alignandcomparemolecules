/***************************************************************************************/ 
/*
 * 
 * Class for output of program  
 * Methods:
 * 	displayDualMatrix -> Display on screen two matrix in a row
 * 							-> Display on screen two vectors in a row
 * 	displayFilesNames -> Display on screen the name of the files which will to use 
 * 	displayItsTheSame -> Display on screen a flag when the main inertia axis are the same in the two molecules 
 * 	saveXYZFile -> Create a file with format xyz for the molecule 
 *
 */
/***************************************************************************************/  
#ifndef _OUTPUT_H_
#define _OUTPUT_H_
#include <string>
using std::string;
#include <vector>
using std::vector;
/***************************************************************************************/ 
#include "screenutils.h"
#include "atomsinmolecule.h"
/***************************************************************************************/ 
/***************************************************************************************/ 

class OutputAlignProgram{
	/***************************************************************************************/ 
	public:
		OutputAlignProgram();

		ScreenUtils scrut;

		void displayDualMatrix(string,vector<vector<double>>,vector<vector<double>>);
		void displayDualMatrix(string,vector<double>,vector<double>);
		void displayFilesNames(bool,string,string);
		void displayResult(bool,string);
		void displayXYZFile(string,vector<Atom>);
		void display_booth_XYZFile(string,string,vector<Atom>,vector<Atom>);
		bool saveXYZFile(string,string,vector<Atom>);

		void displaySameMolecule();
		void displayDifferentMolecule();
		void displayIsomerMolecule();
      void correctInputData(bool,bool);

      void displayRMSD(double);
      void displayLargeEigenVectors(double,vector<vector<double>>,vector<vector<double>>);
      void displayDegeneracy();

	/***************************************************************************************/ 


		void displayInertiaTensorEigenVecEigenVal(vector<vector<double>>,vector<vector<double>> ,vector<vector<double>> ,vector<vector<double>> ,vector<double> ,vector<double>);

	/***************************************************************************************/ 
		// Variables
	/***************************************************************************************/ 
	
	/***************************************************************************************/
	/***************************************************************************************/ 
	private:
		// Variables
	/***************************************************************************************/ 

	/***************************************************************************************/ 
	/***************************************************************************************/ 

	protected:
};
#endif // _OUTPUT_H_
