#ifndef IBETA_
#define IBETA_
/*
 IBETA represents the "beta" track in a bb event (or background)
 eg, a "connected wire" which links to the true electron(s) in the event
 An IBeta is derived from an IBlob (eg a collection of hits) 
 JJGC, July, 2014.

 Notice that we derive from IBlob
 std::vector<const Hit* > GetSortedHits() const {return fSortedHitVector;}

 For an IBeta the GetSortedHits() should give us the collection of hits (sorted in z)
 from which we will fit using kZ
*/


#include "IBlob.h"

// fType values:
// "2e", "1e", "gamma"
namespace alex {

class IBeta: public IBlob {
	public:
		IBeta(){};
		~IBeta();

		void AddBlob(const IBlob& blob);
		
		void SetType(std::string type) {fType = type;}
		std::vector<IBlob*> GetBlobs() const {return fBlobs;}
		std::string GetType() const {return fType;}

	private:
		void ClearBlobs();
		std::vector<IBlob*> fBlobs;
		std::string fType;
			
	};

}
#endif
