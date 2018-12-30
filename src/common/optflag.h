#ifndef _OPTSFLAG_H
#define _OPTSFLAG_H

#include <string>
using std::string;

class OptFlags{
	public:
		OptFlags();
				
		bool display_inertia_tensor;
		bool display_large_eigenvector;
		double vector_increase_length;
		bool display_output_coordenates;
		bool save_output_coordenates;
		bool display_rms;
      bool quiet_version;

		void getOptions(int &argc, char** &argv);
		void printHelpMenu(int &argc, char** &argv); 

};


#endif //_OPTSFLAG_H
