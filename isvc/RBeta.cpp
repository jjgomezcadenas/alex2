
/*
 Hit
 JJGC, July, 2014.
*/

#include "RBeta.h"


namespace alex {

	void RBeta::AddBlob(const IBlob&  blob)
	{
		auto iblob =new IBlob(blob);
		fBlobs.push_back(iblob);
	}

	void RBeta::AddCoreHit(const Hit&  hit)
	{
		auto hitp =new Hit(hit);
		fCoreHit.push_back(hitp);
	}

	void RBeta::AddEffHit(const Hit&  hit)
	{
		auto hitp =new Hit(hit);
		fEffHit.push_back(hitp);
	}

	void RBeta::AddPhoton(const IBeta& phot)
	{
		auto photon =new IBeta(phot);
		fPhotons.push_back(photon);
	}

        void RBeta::ReverseEffHits()
	{
		std::cout << "Reversing: " << fEffHit[0]->XYZ().X() << " and " << fEffHit[fEffHit.size()-1]->XYZ().X() << std::endl;
		std::reverse(fEffHit.begin(),fEffHit.end());
		std::cout << "Now: " << fEffHit[0]->XYZ().X() << " and " << fEffHit[fEffHit.size()-1]->XYZ().X() << std::endl;
	}

	RBeta::~RBeta()
	{
		for(auto iblob : fBlobs)
		  delete iblob; 

		for(auto hit : fCoreHit)
                  delete hit;

                for(auto hit : fEffHit)
                  delete hit;

                for(auto phot : fPhotons)
                  delete phot;

		fCoreHit.clear();
		fEffHit.clear();	
		fBlobs.clear();
		fPhotons.clear();
	}
}
