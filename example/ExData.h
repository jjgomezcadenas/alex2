#ifndef EXDATA_H
#define EXDATA_H

#include <string>
#include <TVector3.h>
#include <alex/IData.h>


//Interface for Alex data. 

namespace alex {

	class ExData: public IData, public INamed {
	public:
		ExData() {};
		virtual ~ExData(){};
		
		std::string Serialize() const ;
		
		void SetData(TVector3 x){fX = x;}
		TVector3 GetData(){return fX;}

		std::string  Name() const {return fName;}
		void SetName(std::string name) {fName = name;}

	private:
	TVector3 fX;
	std::string fName;
	};
}
#endif //IDATA_H