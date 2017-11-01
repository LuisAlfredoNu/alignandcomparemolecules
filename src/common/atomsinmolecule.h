/* 
 * Class for Atoms inside of molecules for align and compare 2 molecules
 */

#include <string>
using std::string;
#include <vector>
using std::vector;
/***************************************************************************************/ 
class Atom{
/***************************************************************************************/ 
public:
	Atom();
/***************************************************************************************/ 
	// Variables
	int atomNumber;
	string atomSymbol;
	vector<double> atomCoordinates; 
	double atomWeight;
	bool statusData;

/***************************************************************************************/ 
/* Assing values to coordinates X,Y,Z */
	void setCoordinates(double x, double y, double z);
/* Get the Value of the atom */
	double getXCoordinate();
	double getYCoordinate();
	double getZCoordinate();
/* Assign name and letter to type of element */
	void setAtomSymbol(string);
	void setAtomNumber(int);

/***************************************************************************************/  
/***************************************************************************************/  

protected:
	double xPosition;
	double yPosition;
	double zPosition;

/***************************************************************************************/  
/***************************************************************************************/  

private:
	string convertAtomNumber2AtomSymbol(int);
	int convertAtomSymbol2AtomNumber(string);
/* Assign atomic weight to the element  */
	double setAtomWeight(int);
};
