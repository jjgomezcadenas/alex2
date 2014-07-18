#ifndef TOYA_H
#define TOYA_H

#include <string>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <alex/IAlgorithm.h>
#include <TH1.h>


//Interface for Alex algorithms. 

namespace alex {

	class ToyAnalysis: public IAlgorithm {
	public:

		bool Init() ;
		bool Execute() ;
		bool End() ;
		std::string  Name(){return fName;}

	private:
		std::string fName;
		TH1F* fH1;
		std::ifstream fIn;

	};
}
#endif 