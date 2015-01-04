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
#include <alex/RBeta.h>
#include "IDefs.h"

#include <TFile.h>
#include <TTree.h>
#include <TRandom2.h>



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
		void CreateTracks(double maxDist);
                void CreateRTracks(double coreDist, int fSparseWidth, double blobRadius);
		void CreateKFObjects(double energy, double errXY, bool addPhotons = false, bool forwardFit = true);
                void GuessInitialMomentum(double energy, std::vector<double> &p0, std::vector<double> xlist, std::vector<double> ylist, std::vector<double> zlist);

		IParticles GetElectrons() const {return fElectrons;}
		int GetNumberOfElectrons() const {return fElectrons.size();}
		IParticles GetPrimaryElectrons() const {return fBetas;}
		int GetNumberOfPrimaryElectrons() const {return fBetas.size();}
		std::pair<IParticle, IParticle> GetPMaxElectrons() const {return fBetasMax;}
		IHits GetTrueHits() const {return fTrueHits;}
		std::pair<IHits, IHits> GetPMaxElectronsHits() const {return fBetasMaxHits;}
		std::vector<const IBeta*> GetIBetas() const { return fIBetas; }
                RBeta* GetRBeta() const { return fRBeta; }
	        std::vector<double> GetMinDistList() const { return minDistList; }

                std::vector<const Hit*> GetKFHits() const { return fKFHits; }
                std::vector<double> GetKFMErrors() const { return fKFMErrors; }
                std::vector<double> GetKFv0() const { return fKFv0; }
                std::vector<double> GetKFp0() const { return fKFp0; }

                void SetStartEvent(int sevt) { startEvt = sevt; }
                int GetStartEvent() const { return startEvt; }

                void SetEvtNum(int sevt) { evtNum = sevt; }
                int GetEvtNum() const { return evtNum; }

	private:
		void FetchElectrons();
		void FetchPMaxElectrons();
                double ComputeMinDist(IHits ihs1, IHits ihs2);
		double ComputeMinDist(std::vector<const alex::Hit*> ihs1, IHit ih, int& imin, int& iside);
                std::vector<double> ApplyLPF(const std::vector<double> & bc, const std::vector<double> & ac, const std::vector<double> & xv);

        int startEvt;
        int evtNum;

        TRandom2* fRandom;
	TFile* fIfile;
  	TTree* fEvtTree;
        RBeta* fRBeta;
  	const irene::Event* fIevt;
  	IParticles fElectrons;
  	IParticles fBetas; //beta = primary electron
  	IHits fTrueHits;
  	std::vector<const irene::Track*> fIreneTracks ; 
  	std::pair<IParticle, IParticle> fBetasMax;
  	std::pair<IHits, IHits> fBetasMaxHits;
        std::vector<const IBeta*> fIBetas;
        std::vector<double> minDistList;
        
        std::vector<const Hit*> fKFHits;
        std::vector<double> fKFMErrors;
        std::vector<double> fKFv0;
        std::vector<double> fKFp0;
			
	};

	typedef alex::SingletonTemplate<IreneManager> ISvc;   // Global declaration

}
#endif
