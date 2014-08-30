#ifndef HITXX_
#define HITXX_
/*
 A hit
*/

#include "IDefs.h"
#include <TVector3.h>

namespace alex {

class Hit {
	public:
		Hit(const IHit&  ihit);
		Hit(const Hit& hit);
		~Hit(){};

		double Edep() const {return fEne;}
		TVector3 XYZ()const {return fXYZ;}
		std::string Type() const {return fType;}
		int Owner() const {return fOwner;}
		void SetType(std::string type){fType=type;}
		void SetOwner(int owner){fOwner=owner;}


	private:
	std::string fType; // true, sparsed, smeared
	TVector3 fXYZ; 
	double fEne; // (x,y,z,e)
	int fOwner; // (1)= beta 1; (2) = beta 2; (0) not beta 1 or beta 2
	};

}
#endif