
/*
 Hit
 JJGC, July, 2014.
*/

#include "IBlob.h"


namespace alex {

	IBlob::IBlob(const IBlob&  iblob)
	{
		std::vector<const Hit* > ihits =iblob.GetHit();		
		for(auto hit : ihits)
			this->AddHit(*hit);
	}

	void IBlob::AddHit(const Hit&  hit)
	{
		auto hitp =new Hit(hit);
		fHitVector.push_back(hitp);
	}

	void IBlob::AddHit(const IHit& hit)
	{
		auto hitp = new Hit(hit);
		fHitVector.push_back(hitp);
	}

  void IBlob::Clear()
  {
    fHitVector.clear();
  }

  void IBlob::DeleteHits()
  {
    for(auto hit : fHitVector)
    	delete hit; 
		fHitVector.clear();	
  }

	IBlob::~IBlob()
	{
		DeleteHits();	
	}
}
