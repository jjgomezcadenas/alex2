#ifndef IBLOB_
#define IBLOB_
/*
 IBLOB represents a blob at the end of the IBeta track.
 A blob is a set of hits 
 JJGC, JR, October, 2014.
*/


#include "Hit.h"

namespace alex {

class IBlob {
	public:
		IBlob(){};
		IBlob(const IBlob&  iblob);
		virtual ~IBlob();
		void AddHit(const Hit& hit);
		void AddHit(const IHit& hit);
		
		std::vector<const Hit* > GetHit() const {return fHitVector;}

	private:
		void Clear();
		void DeleteHits();
		std::vector<const Hit*> fHitVector;
			
	};

}
#endif
