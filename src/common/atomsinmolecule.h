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
/* Assing values to coordinates X,Y,Z */
	void setCoordinates(double x, double y, double z);
/* Get the Value of the atom */
	double getXCoordinate();
	double getYCoordinate();
	double getZCoordinate();
	//vector<double> getCoordinates();
/* Assign name and letter to type of element */
	void setTypeElement(string element);
	string getTypeElement();
/* Assign atomic weight to the element  */
	double getWeightElement();

/***************************************************************************************/  
/***************************************************************************************/  

protected:
	double xPosition;
	double yPosition;
	double zPosition;
	vector<double> coordinates; 
	string typeAtom;
	int numelement;
	double weighelement;

/***************************************************************************************/  
/***************************************************************************************/  

private:
	int convertTypeElement2NumElement(string);
};
