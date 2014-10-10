
/*
 Hit
 JJGC, July, 2014.
*/

#include "IBeta.h"


namespace alex {


	void IBeta::AddHit(const Hit&  hit)
	{
		auto hitp =new Hit(hit);
		fHitVector.push_back(hitp);
	}

	void IBeta::AddHit(const IHit& hit)
	{
		auto hitp = new Hit(hit);
		fHitVector.push_back(hitp);
	}

        void IBeta::Clear()
        {
                fHitVector.clear();
        }

	IBeta::~IBeta()
	{
		for(auto hit : fHitVector)
    {
      delete hit; 
		}	
		fHitVector.clear();	
	}
}
