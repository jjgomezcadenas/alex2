#ifndef TOYDATA_H
#define TOYDATA_H

#include <string>
#include <TVector3.h>
#include <alex/IData.h>


//Interface for Alex data. 

namespace alex {

	class ToyData: public IData {
	public:
		ToyData() {};
		virtual ~ToyData(){};
		
		std::string Serialize() const ;
		void Recreate(std::string) ;
	
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