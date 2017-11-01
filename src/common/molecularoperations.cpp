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
#endif // _MOLECULAR_OPERATIONS_CPP_

