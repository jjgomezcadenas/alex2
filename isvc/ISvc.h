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

#include <TFile.h>
#include <TTree.h>

#include <irene/Event.h>
#include <irene/Track.h>
#include <irene/Particle.h>
#include <TLorentzVector.h>


namespace alex {
typedef std::pair<TLorentzVector,double> IHit;
typedef std::vector<std::pair<TLorentzVector,double> > IHits;
typedef const irene::Particle* IParticle;
typedef std::vector<const irene::Particle*> IParticles;

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
  	std::vector<const irene::Track*> fIreneTracks ; 
			
	};

	typedef SingletonTemplate<IreneManager> ISvc;   // Global declaration

}
#endif