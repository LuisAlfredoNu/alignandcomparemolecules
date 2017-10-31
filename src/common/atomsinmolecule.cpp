/* 
 * Class for Atoms inside of molecules for align and compare 2 molecules
 */

#include <iostream>
using std::cout;
using std::endl;
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
	numelement=0;
	weighelement=0.0;
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
vector<double> Atom::getCoodinates(){
	return coordinates
}
/***************************************************************************************/ 
/***************************************************************************************/ 
void Atom::setTypeElement(string element){
	if(element.size()<3){
		typeAtom = element;}
	else{
		typeAtom = element; 
		typrAtom +=" Invalid element";}
			// Launch a error menssage
}
string Atom::getTypeElement(){
	return typeAtom;
}
/***************************************************************************************/  
/***************************************************************************************/  
double Atom::getWeightElement(){
	switch(convertTypeElement2NumElement(typeAtom)){
		case 1 :
			weighelement=1.002;
			break;
		case 0 :
			weighelement=0.0;
			// Launch a error menssage
			break;
	}
	return weighelement;
}
int Atom::convertTypeElement2NumElement(string element){
	if(element == "H"){
		numelement=1;}
	else{numelement = 0;}
			// Launch a error menssage
	return numelement;
}
