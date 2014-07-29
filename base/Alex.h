#ifndef ACONFIG_
#define ACONFIG_
/*
 ACONFIG: Alex Configuration
 

 JJGC, July, 2014.
*/

#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <map>
#include <alex/SingletonTemplate.h>
#include <alex/IAlgorithm.h>
 #include <alex/IData.h>

#include <alex/AConfigSvc.h>
#include <alex/LogUtil.h>

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>


namespace alex {

class AlexManager {
	public:
		AlexManager(){};
		virtual ~AlexManager(){};
		void Init(std::string debugLevel);

		void RegisterAlgorithm(IAlgorithm* algo );
		void InitAlgorithms();
		void ExecuteAlgorithms();
		void EndAlgorithms();
		void ClearAlgorithms();

		void RegisterData(std::string name, IData* data )
		{fIData[name]=data;}
		IData* RetrieveData(std::string name){return fIData[name]; }
		void ClearData();
		// const IData* RetrieveData(std::string name) const
		// {return fIData[name];} 

		void InitHistoFile(std::string fileName);
		void WriteHistoFile();
		void CloseHistoFile();

	private:
		std::vector<IAlgorithm*> fIAlgo;
		std::map<std::string,IData*> fIData;
		
		TFile* fHistoFile;
			
	};

	typedef SingletonTemplate<AlexManager> Alex;   // Global declaration

}
#endif