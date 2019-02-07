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
	display_rms = false;
   quiet_version = false;
   run_big_test = false;
};
/***************************************************************************************/ 
void OptFlags::getOptions(int &argc, char** &argv){
	ScreenUtils scrut;
	if (argc<2) {
      cout << endl;
      scrut.PrintScrStarLine();
      scrut.DisplayErrorMessage("Two input files are needed!!! ");
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
				case 'r':
					display_rms = true;
					break;
				case 'B':
					run_big_test = true;
					break;
            case 'q':
               quiet_version = true;
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
   if(quiet_version){
      display_inertia_tensor = false;
      display_large_eigenvector = false;
      display_output_coordenates = false;
      save_output_coordenates = false;
      display_rms = false;
      run_big_test = false;
   }
}
/***************************************************************************************/  
void OptFlags::printHelpMenu(int &argc, char** &argv){
   string progname=argv[0];
   size_t pos=progname.find("./");
   if (pos!=string::npos) {progname.erase(pos,2);}

	cout << endl;
   cout << "\nUsage:\n\n\t" << progname << " -flag ... -flag moleculeA.xyz moleculeB.xyz \n\n";
   cout << "Where moleculeA.xyz and moleculeB.xyz contain the input structures. The files must be of type XYZ, and options can be:\n\n"
        << "  -m        \tDisplays the inertia tensor on screen, the eigenvectors of the inertia tensor and its eigenvalues." << endl
        << "  -c        \tDisplays the coordinates of molecules after all alignment processes have been performed. " << endl
        << "  -s        \tSaves the coordinates of the atoms at the end of the processes in XYZ file with \"_align.xyz\" appended to the file name. " << endl
        << "  -r        \tDisplays the RMSD between molecules" << endl
        << "  -q        \tQuiet output version. Displays the result, without color and in a single line." << endl;
   cout << "  -h        \tDisplay the help menu.\n\n";
}




#endif //_OPTSFLAG_CPP
