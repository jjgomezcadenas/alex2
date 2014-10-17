#ifndef RBETA_
#define RBETA_
/*
 RBETA represents the "reconstructed" beta... it is an IBeta with a log of s... sense

 Notice that we derive from IBeta
 
*/


#include "IBeta.h"

namespace alex {

class RBeta: public IBeta {
	public:
		RBeta(){};
		virtual ~RBeta();

		void AddBlob(const IBlob& blob);
		void AddCoreHit(const Hit& hit);
		void AddEffHit(const Hit& hit);
		void AddPhoton(const IBeta& photon);

                void ReverseEffHits();
		
		std::vector<const IBlob*> GetBlobs() const {return fBlobs;}
		std::vector<const IBeta*> GetPhotons() const {return fPhotons;}
		std::vector<const Hit*> GetCoreHits() const {return fCoreHit;}
		std::vector<const Hit*> GetEffHits() const {return fEffHit;}
		
	private:
		std::vector<const Hit*> fCoreHit;
		std::vector<const Hit*> fEffHit;
		std::vector<const IBlob*> fBlobs;
		std::vector<const IBeta*> fPhotons;
		
			
	};

}
#endif
