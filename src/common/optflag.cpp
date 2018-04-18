#ifndef _OPTSFLAG_CPP
#define _OPTSFLAG_CPP

#include "optflag.h"
#include "screenutils.h"

#include <iostream>
using std::cout;
using std::endl;
using std::ios;
#include <cstdlib>
using std::exit;
#include <string>
using std::string;
/***************************************************************************************/  
OptFlags::OptFlags()	{
	display_inertia_tensor = false;
	display_large_eigenvector = false;
	vector_increase_length = 0.0;
	display_output_coordenates =false;
	save_output_coordenates = false;
};
/***************************************************************************************/ 
void OptFlags::getOptions(int &argc, char** &argv){
	ScreenUtils scrut;
	if (argc<2) {
      cout << endl;
      scrut.PrintScrStarLine();
      scrut.DisplayErrorMessage("Dont input files in arguments ");
      cout << "\t Example:" << argv[0] <<" file_A.xyz file_B.xyz  "<< endl;
      scrut.PrintScrStarLine();
      cout << endl;
		exit(1);
	}
	for (int i=1; i<argc; i++){
		if (argv[i][0] == '-'){
			switch (argv[i][1]){
				case 'm':
					display_inertia_tensor = true;
					break;
				case 'l':
					display_large_eigenvector = true;
					sscanf(argv[++i],"%lf",&vector_increase_length);
					break;
				case 'c':
					display_output_coordenates = true;
					break;
				case 's':
					save_output_coordenates = true;
					break;
				case 'h':
					printHelpMenu(argc,argv);
					exit(1);
					break;
				default:
					cout << "\nCommand line error. Unknown switch: " << argv[i] << endl;
					cout << "\nTry: \n\t" << argv[0] << " -h\n" << endl << "to view the help menu.\n\n";
					exit(1);
			}
		}
	}
}
/***************************************************************************************/  
void OptFlags::printHelpMenu(int &argc, char** &argv){
   string progname=argv[0];
   size_t pos=progname.find("./");
   if (pos!=string::npos) {progname.erase(pos,2);}

	cout << endl;
   cout << "\nUsage:\n\n\t" << progname << " -flag ... -flag moleculeA.xyz moleculeB.xyz \n\n";
   cout << "Where moleculeA.xyz and moleculeB.xyz is the input structure in file type XYZ, \n\t and options can be:\n\n"
        << "  -m        \tDisplay on screen the inertia tensor, eigenvectors of inertia tensor" << endl
		  << "            \t  and eigenvalues of inertia tensor." << endl
        << "  -l n      \tDisplay on screen the coordinates of eigenvector with n increase of them length " << endl
        << "  -c        \tDisplay on screen the coordinates of atoms after all procces " << endl
        << "  -s        \tSave the coordinates of atoms after all procces in XYZ file with *_align.xyz terminataion " << endl;
   cout << "  -h        \tDisplay the help menu.\n\n";
}




#endif //_OPTSFLAG_CPP
