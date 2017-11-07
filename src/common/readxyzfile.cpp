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
#include <string>
using std::string;
/***************************************************************************************/ 
#include "readxyzfile.h"
/***************************************************************************************/  
/***************************************************************************************/  

ReadXYZFile::ReadXYZFile(){ 
	open_without_problems = true;
}

bool ReadXYZFile::getValuesFromFile(string filename, vector<Atom> & molecule){

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

	molecule.resize(getNumofAtoms(XYZFile));
	getDataAtoms(XYZFile,molecule);


	XYZFile.close();

	return open_without_problems;
}
/***************************************************************************************/ 
int ReadXYZFile::getNumofAtoms(ifstream &file){

	int numof_Atoms=0;

	file.seekg(0,file.beg);

	file >> numof_Atoms;

	if(numof_Atoms == 0) open_without_problems = false;

	return numof_Atoms;
}
/***************************************************************************************/ 
bool ReadXYZFile::getDataAtoms(ifstream &file, vector<Atom>& molecule){

	int data_startpos;
	int i = 0;
	int numatoms = getNumofAtoms(file);
	
	int atomnumber;
	string symbol;
	vector<double> vectorposition (3);
	
	if(typeDataNumOChar(file)){
		while(i < numatoms ){
			file.seekg(begindata_pos);
			file >> atomnumber;
			file >> vectorposition[0];
			file >> vectorposition[1];
			file >> vectorposition[2];

			molecule[i].setAtomNumber(atomnumber);
			molecule[i].setCoordinates(vectorposition);

			cout << "count  i = " << i << "begindata_pos" << begindata_pos << endl;
			i++;
		}
	}else{
		while(i < numatoms){
			file.seekg(0);
			getline(file,symbol);
		//	file >> symbol;
		//	file >> vectorposition[0];
		//	file >> vectorposition[1];
		//	file >> vectorposition[2];

			cout << "Symbol" <<symbol << endl;
		//	molecule[i].setAtomSymbol(symbol);
		//	molecule[i].setCoordinates(vectorposition);
			
			cout << "count  i = " << i << endl;
			cout << "count  i = " << i << "begindata_pos" << begindata_pos << endl;

			i++;
		}
	}

}
/***************************************************************************************/ 
bool ReadXYZFile::typeDataNumOChar(ifstream &file){

	bool is_number = true;
	int number;
	string line;
	
	file.seekg(0,file.beg);

	for(int i=0; i<2;i++) getline(file,line);

	file >> number;
	begindata_pos = file.tellg();
			cout << "begindata_pos" << begindata_pos << endl;

	if(file.fail()) is_number=false;
	
	return is_number;
}
/***************************************************************************************/ 
/***************************************************************************************/ 
#endif // _READXYZFILE_CPP_

