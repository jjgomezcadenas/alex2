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


namespace alex {

class IreneManager {
	public:
		IreneManager(){};
		virtual ~IreneManager(){};
		void Init(std::string debugLevel);
		void InitDstFile(std::string fileName);
		int DstEntries();
		int DstGetEntry(int ivt);
		void LoadEvent(const irene::Event* ievt);
		const irene::Event& GetEvent();

	private:
		
		TFile* fIfile;
  	TTree* fEvtTree ;
  	const irene::Event* fIevt;
			
	};

	typedef SingletonTemplate<IreneManager> ISvc;   // Global declaration

}
#endif