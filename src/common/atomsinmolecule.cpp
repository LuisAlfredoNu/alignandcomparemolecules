/* 
 * Class for Atoms inside of molecules for align and compare 2 molecules
 */

#include <string>
using std::string;
#include <vector>
using std::vector;
#include "atomsinmolecule.h"
/***************************************************************************************/ 
/***************************************************************************************/ 

Atom::Atom(){
	coordinates.assign(3,0.0);
	xPosition=0.0;
	yPosition=0.0;
	zPosition=0.0;
	typeAtom.clear();
	}

void Atom::setCoordinates(double x, double y, double z){
	xPosition=x;
	yPosition=y;
	zPosition=z;
	//coordinates[0]=x;
	//coordinates[1]=y;
	//coordinates[2]=z;
}
/***************************************************************************************/ 
double Atom::getXCoordinate(){
	return xPosition;
}
double Atom::getYCoordinate(){
	return yPosition;
}
double Atom::getZCoordinate(){
	return zPosition;
}
//vector<double> Atom::getCoodinates(){
//	return coordinates
//}
/***************************************************************************************/ 
/***************************************************************************************/ 
void Atom::setTypeElement(string element){
	if(element.size()>3){
		typeAtom = element;}
	else{
		cout << "Its not a element" << endl;
}
string Atom::getTypeElement(){
	return typeAtom;
}
