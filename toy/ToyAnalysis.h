#ifndef TOYA_H
#define TOYA_H


#include <alex/IAlgorithm.h>
#include <TH1.h>


//Interface for Alex algorithms. 

namespace alex {

	class ToyAnalysis: public IAlgorithm {
	public:

		bool Init() ;
		bool Execute() ;
		bool End() ;
		
		std::string  Name() const {return fName;}
		void SetName(std::string name) {fName = name;}


	private:
		std::string fName;
		TH1F* fH1;

	};
}
#endif 