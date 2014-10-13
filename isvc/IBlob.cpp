
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
		fBlob = new Hit(*(iblob.GetBlob()));
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

	void IBlob::CreateBlob(double x, double y, double z, double E)
	{
		fBlob = new Hit(x, y, z, E);

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
		delete fBlob;
	}
}
