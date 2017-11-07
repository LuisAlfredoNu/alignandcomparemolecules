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
	Natoms = 0;
	begindata_pos = 0;
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

	Natoms = getNumofAtoms(XYZFile);
	molecule.resize(Natoms);
	getDataAtoms(XYZFile,molecule);
	open_without_problems = statusAllData(molecule);

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
void ReadXYZFile::getDataAtoms(ifstream &file, vector<Atom>& molecule){

	int atomnumber;
	string symbol;
	vector<double> vectorposition (3);
			
	file.seekg(begindata_pos);
	
	int i = 0;
	if(typeDataNumOChar(file)){
		while(i < Natoms ){
			file >> atomnumber;
			file >> vectorposition[0];
			file >> vectorposition[1];
			file >> vectorposition[2];

			molecule[i].setAtomNumber(atomnumber);
			molecule[i].setCoordinates(vectorposition);

			i++;
		}
	}else{
		while(i < Natoms){
			file >> symbol;
			file >> vectorposition[0];
			file >> vectorposition[1];
			file >> vectorposition[2];

			molecule[i].setAtomSymbol(symbol);
			molecule[i].setCoordinates(vectorposition);
			
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
	
	begindata_pos = file.tellg();		
	cout << "begindata_pos" << begindata_pos << endl;

	file >> number;
	if(file.fail()){
		is_number = false;
		file.clear();
	}
	
	return is_number;
}
/***************************************************************************************/ 
bool ReadXYZFile::statusAllData(vector<Atom> molecule){

	for(int i=0;i<Natoms;i++)
		if(!molecule[i].statusData) return false;

	return true;
}
/***************************************************************************************/ 
/***************************************************************************************/ 
#endif // _READXYZFILE_CPP_

