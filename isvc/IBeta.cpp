
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

	IBeta::~IBeta()
	{
		for(auto hit : fHitVector)
    {
      delete hit; 
		}	
		fHitVector.clear();	
	}
}