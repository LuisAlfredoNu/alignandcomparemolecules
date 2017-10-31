/***************************************************************************************/  
/* Class of Molecular Operations       */
/***************************************************************************************/  

#include "atomsinmolecule.h"
/***************************************************************************************/  
/***************************************************************************************/  

MolecularOpetations::MolecularOpetations(){ }


vector<double> MolecularOpetations::massCenter(vector<Atom> molecule){
	double averagemass=0.0;

	for(unsigned int i=0;i< molecule.size(); ++i){
		averagemass += molecule[i].getWeightElement();
	}
	vector<double> rmasscenter (3,0.0);
	for(unsigned int i=0;i<3;++i){
		for(unsigned int ii=0;i< molecule.size();++i){
			rmasscenter[i]=molecule[ii].getWeightElement()*molecule[ii].getXCoordinate();

		}
	}
}
