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
		Hit(double x, double y, double z, double E);
                Hit(double x, double y, double z, double zend, double E);
		~Hit(){};

		double Edep() const {return fEne;}
                double MinDist() const {return fMinDist;}
		TVector3 XYZ()const {return fXYZ;}
		std::string Type() const {return fType;}
		int Owner() const {return fOwner;}
		void SetType(std::string type){fType=type;}
		void SetOwner(int owner){fOwner=owner;}
                void SetMinDist(double minDist) {fMinDist = minDist;}
		void AddEnergy(double ene) { fEne += ene; }

	private:
                std::string fType; // true, sparsed, smeared
                TVector3 fXYZ;
                double zend; // end z for Hit representing a slice
                double fEne; // (x,y,z,e)
                int fOwner; // (1)= beta 1; (2) = beta 2; (0) not beta 1 or beta 2
                double fMinDist; // holds minimum distance to main track, if applicable
	};

}
#endif
