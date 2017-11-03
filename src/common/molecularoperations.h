
/***************************************************************************************/  
/* Class of Molecular Operations
 *	Class for making a different operations to arrays of atoms -> molecule 
 *		massCenter -> Find the center of mass of the moelcules
 *		inertiaTensor -> Make the inertia tensor of the molecule 
 *		moveCM2Origin -> Move the molecule where the center of mass and origin of coordinate system are the same 
 *
 * */
/***************************************************************************************/  
/***************************************************************************************/  
#ifndef _MOLECULAR_OPERATIONS_H_
#define _MOLECULAR_OPERATIONS_H_

#include <vector>
using std::vector;
/***************************************************************************************/ 
#include "atomsinmolecule.h"

class MolecularOperations{
	/***************************************************************************************/ 
	public:
		MolecularOperations();
	/***************************************************************************************/ 
		// Variables
		vector<double> centerMass;
	/***************************************************************************************/
		vector<double> massCenter(vector<Atom>);
		vector<vector<double>> inertiaTensor(vector<Atom>);
		vector<Atom> moveCM2Origin(vector<Atom>);
	/***************************************************************************************/ 
	/***************************************************************************************/ 
	private:

	/***************************************************************************************/ 
	/***************************************************************************************/ 

	protected:
};
#endif // _MOLECULAR_OPERATIONS_H
