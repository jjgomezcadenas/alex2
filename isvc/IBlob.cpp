
/*
 Hit
 JJGC, July, 2014.
*/

#include "IBlob.h"


namespace alex {

	IBlob::IBlob()
	{
		fBlob = new Hit(0.,0.,0.,0.);
	}

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

	void IBlob::InsertHit(const Hit& hit, int i)
	{
		auto hitp = new Hit(hit);
		fHitVector.insert(fHitVector.begin() + i, hitp);
	}

	void IBlob::InsertHit(const IHit& hit, int i)
	{
		auto hitp = new Hit(hit);
		fHitVector.insert(fHitVector.begin() + i, hitp);
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

  double IBlob::GetEnergy() const
  {
    double e = 0.;
    for(auto hit : fHitVector)
      e += hit->Edep();
    return e;
  }

	IBlob::~IBlob()
	{
		DeleteHits();	
		delete fBlob;
	}
}
