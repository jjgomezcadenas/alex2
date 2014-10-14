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
// "LT = long track (> 200 hits) ", "ST = (> 50 < 200) ", "PH= photon < 50"
// fMCType: --> "1e", "2e "gamma"
namespace alex {

class IBeta: public IBlob {
	public:
		IBeta(){};
		virtual ~IBeta(){};
		IBeta(const IBeta& ibeta);

		void SetType(std::string type) {fType = type;}
		std::string GetType() const {return fType;}

		void SetMCType(std::string type) {fType = type;}
		std::string GetMCType() const {return fMcType;}

		void SetV0(TVector3 V0) {fV0 = V0;}
		TVector3 GetV0() const {return fV0;}

		void SetVX(TVector3 V1,TVector3 V2);
		std::pair<TVector3, TVector3> GetVX() const {return fVx;}

	private:
		std::string fType;
		std::string fMcType;
		TVector3 fV0; // the start (vertex of the 1e or 2e)
		std::pair<TVector3, TVector3> fVx; // for 1 e fVx.first is the end
		// for 2e fVx.first and fVx.second are the two ends. 
			
	};

}
#endif
