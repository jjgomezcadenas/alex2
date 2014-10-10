#ifndef IBETA_
#define IBETA_
/*
 IBETA represents the "beta" track in a bb event (or background)
 eg, a "connected wire" which links to the true electron(s) in the event 
 JJGC, July, 2014.
*/


#include "Hit.h"
#
namespace alex {

class IBeta {
	public:
		IBeta(){};
		~IBeta();
		void AddHit(const Hit& hit);
		void AddHit(const IHit& hit);
		void Clear();
		void SetType(std::string type) {fType = type;}
		std::vector<const Hit* > GetHit() const {return fHitVector;}
		std::string GetType() const {return fType;}

	private:
		std::vector<const Hit*> fHitVector;
		std::string fType;
			
	};

}
#endif
