/***************************************************************************************/  
/* Class of ReadXYZFile       */
/***************************************************************************************/  
#ifndef _READXYZFILE_CPP_
#define _READXYZFILE_CPP_

#include <iostream>
using std::cout;
using std::endl;
using std::ios;
#include <vector>
using std::vector;
/***************************************************************************************/ 
#include "readxyzfile.h"
/***************************************************************************************/  
/***************************************************************************************/  

ReadXYZFile::ReadXYZFile(){ 
	open_without_problems = true;
}

bool ReadXYZFile::getValuesFromFile(string filename, vector<Atom> molecule){

	ifstream XYZFile;

	XYZFile.open(filename.c_str(),ios::in);

	if (!(XYZFile.good())) {
		cout << "Error: File " << filename << "could not be opened...\n";
#if DEBUG
		cout << __FILE__ << ", line: " << __LINE__ << endl;
#endif
		open_without_problems = false;
		return false;
	}

	molecule.resize(getNumofAtoms());

	return open_without_problems;
}
/***************************************************************************************/ 
int getNumofAtoms(ifstream file){

	int numof_Atoms=0;

	return numof_Atoms;
}

/***************************************************************************************/ 

/***************************************************************************************/ 
/***************************************************************************************/ 
#endif // _READXYZFILE_CPP_

