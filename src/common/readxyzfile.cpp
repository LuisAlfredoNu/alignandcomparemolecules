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
	sortingAtoms(molecule);
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

	int atomnumber=0;
	string symbol;
	vector<double> vectorposition (3);
			
	
	int i = 0;
	if(typeDataNumOChar(file)){
		file.seekg(begindata_pos);
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
		file.seekg(begindata_pos);
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
void ReadXYZFile::sortingAtoms(vector<Atom>& molecule){

	// Here I adapt the Comb sort algorithm from Wikipedia  
	Atom atom4swap_tmp;

	int swaps = 1;
	int gap = molecule.size();
	int sizemolecule = molecule.size();

	while(!(swaps == 0 && gap ==1)){
		if(gap > 1){
			gap /= 1.3;
			if (gap ==10 || gap == 9){
				gap = 11;
			}
		}

		int i=0; swaps=0;
		while(i + gap < sizemolecule){
			//Sorting by X
			if(molecule[i].atomCoordinates[2] > molecule[i+gap].atomCoordinates[2]){

				atom4swap_tmp = molecule[i];
				molecule[i] = molecule[i+gap];
				molecule[i+gap] = atom4swap_tmp;

				swaps += 1; 
			}else{
				//Sorting by Y
				if(molecule[i].atomCoordinates[2] == molecule[i+gap].atomCoordinates[2]){
					if(molecule[i].atomCoordinates[1] > molecule[i+gap].atomCoordinates[1]){

						atom4swap_tmp = molecule[i];
						molecule[i] = molecule[i+gap];
						molecule[i+gap] = atom4swap_tmp;

						swaps += 1; 
					}else{
						// Sorting by Z
						if(molecule[i].atomCoordinates[1] == molecule[i+gap].atomCoordinates[1]){
							if(molecule[i].atomCoordinates[0] > molecule[i+gap].atomCoordinates[0]){

								atom4swap_tmp = molecule[i];
								molecule[i] = molecule[i+gap];
								molecule[i+gap] = atom4swap_tmp;

								swaps += 1; 
							}
						}
					}
				}
			}
			i += 1;
		}
	}
}
/***************************************************************************************/ 
/***************************************************************************************/ 
#endif // _READXYZFILE_CPP_

