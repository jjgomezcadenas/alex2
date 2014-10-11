#ifndef ISVC_
#define ISVC_
/*
 Irene Svc to provide access to irene (eg DST) to algos
 JJGC, July, 2014.
*/

#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <map>
#include <alex/SingletonTemplate.h>
#include <alex/LogUtil.h>
#include <alex/IBeta.h>
#include "IDefs.h"

#include <TFile.h>
#include <TTree.h>




namespace alex {

class IreneManager {
	public:
		IreneManager(){};
		virtual ~IreneManager(){};
		void Init(std::string debugLevel);
		void InitDst(std::string fileName,const irene::Event* ievt);
		int DstEntries();
		int DstGetEntry(int ivt);
		void LoadEvent(const irene::Event* ievt);
		const irene::Event& GetEvent();
		void CreateTracks();
                //int GetIMaxIBeta();

		IParticles GetElectrons() const {return fElectrons;}
		int GetNumberOfElectrons() const {return fElectrons.size();}
		IParticles GetPrimaryElectrons() const {return fBetas;}
		int GetNumberOfPrimaryElectrons() const {return fBetas.size();}
		std::pair<IParticle, IParticle> GetPMaxElectrons() const {return fBetasMax;}
		IHits GetTrueHits() const {return fTrueHits;}
		std::pair<IHits, IHits> GetPMaxElectronsHits() const {return fBetasMaxHits;}
		std::vector<const IBeta*> GetIBetas() const { return fIBetas; }
        std::vector<double> GetMinDistList() const { return minDistList; }

	private:
		void FetchElectrons();
		void FetchPMaxElectrons();
                double ComputeMinDist(IHits ihs1, IHits ihs2);
		double ComputeMinDist(std::vector<const alex::Hit*> ihs1, IHit ih);

	double fMaxDist;
	TFile* fIfile;
  	TTree* fEvtTree ;
  	const irene::Event* fIevt;
  	IParticles fElectrons;
  	IParticles fBetas; //beta = primary electron
  	IHits fTrueHits;
  	std::vector<const irene::Track*> fIreneTracks ; 
  	std::pair<IParticle, IParticle> fBetasMax;
  	std::pair<IHits, IHits> fBetasMaxHits;
        std::vector<const IBeta*> fIBetas;
        std::vector<double> minDistList;
			
	};

	typedef SingletonTemplate<IreneManager> ISvc;   // Global declaration

}
#endif
