/***************************************************************************************/  
/* Class of Molecular Operations       */
/***************************************************************************************/  
#ifndef _MOLECULAR_OPERATIONS_CPP_
#define _MOLECULAR_OPERATIONS_CPP_

#include <vector>
using std::vector;
/***************************************************************************************/ 
#include "atomsinmolecule.h"
#include "molecularoperations.h"
/***************************************************************************************/  
/***************************************************************************************/  

MolecularOperations::MolecularOperations(){ }


vector<double> MolecularOperations::massCenter(vector<Atom> molecule){
	double totalmass=0.0;

	for(unsigned int i=0;i< molecule.size(); ++i){
		totalmass += molecule[i].atomWeight;
	}
	vector<double> rmasscenter (3,0.0);

	for(unsigned int i=0;i<3;++i){
		for(unsigned int ii=0;ii< molecule.size();++ii){
			rmasscenter[i] += molecule[ii].atomWeight * molecule[ii].atomCoordinates[i];
		}
		rmasscenter[i] /= totalmass;
	}
	return rmasscenter;
}
vector<vector<double>> MolecularOperations::inertiaTensor(vector<Atom> molecule){

	vector<vector<double>> matrix(3,vector<double>(3,0));
	double module_r2=0.0;

	for(int i=0;i<3;i++){
		for(int j=0;j<3;j++){
			if(i==j){
				for(unsigned int im=0;im<molecule.size();im++){
					module_r2 = 0.0;
					for(int xyz=0;xyz<3;xyz++) module_r2 += molecule[im].atomCoordinates[xyz] * molecule[im].atomCoordinates[xyz];
					matrix[i][j] += molecule[im].atomWeight * (module_r2 - molecule[im].atomCoordinates[i] * molecule[im].atomCoordinates[j]); 
				}
			}else{
				for(unsigned int im=0;im<molecule.size();im++) 
					matrix[i][j] += molecule[im].atomWeight * molecule[im].atomCoordinates[i] * molecule[im].atomCoordinates[j] ;
				matrix[i][j] = -matrix[i][j];
			}
		}
	}
	return matrix;
}
/***************************************************************************************/ 
/***************************************************************************************/ 
#endif // _MOLECULAR_OPERATIONS_CPP_

