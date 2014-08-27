#ifndef IBETA_
#define IBETA_
/*
 IBETA represents the "beta" track in a bb event (or background)
 eg, a "connected wire" which links to the true electron(s) in the event 
 JJGC, July, 2014.
*/


#include "Hit.h"
#
namespace alex {

class IBeta {
	public:
		IBeta(){};
		~IBeta(){};
		void AddHit(const Hit& hit);
		void InitDst(std::string fileName,const irene::Event* ievt);
		int DstEntries();
		int DstGetEntry(int ivt);
		void LoadEvent(const irene::Event* ievt);
		const irene::Event& GetEvent();

		IParticles GetElectrons() const {return fElectrons;}
		int GetNumberOfElectrons() const {return fElectrons.size();}
		IParticles GetPrimaryElectrons() const {return fBetas;}
		int GetNumberOfPrimaryElectrons() const {return fBetas.size();}
		std::pair<IParticle, IParticle> GetPMaxElectrons() ;

	private:
		void FetchElectrons();
		
		TFile* fIfile;
  	TTree* fEvtTree ;
  	const irene::Event* fIevt;
  	IParticles fElectrons;
  	IParticles fBetas; //beta = primary electron
  	IHits fTrueHits;
			
	};

	typedef SingletonTemplate<IreneManager> ISvc;   // Global declaration

}
#endif