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
#include <tinyxml2.h>


namespace alex {

class AlexConf {
	public:
		AlexConf(){};
		virtual ~AlexConf(){};
		void Init(std::string debugLevel, std::string rootName);
		void ParseConfiguration(std::string configFile);
		
		std::string SerializeAConfHeader() const;
		std::string SerializeAConfCPP() const;
		std::string SerializeAlgoNames() const;
		std::string SerializeAlgoPaths() const;

	private:
		void ParseAlgosConfiguration();
		void ParseAlgos();
		//tinyxml2::XMLElement* ParseRoot() ;
		std::pair<std::string,std::string> ParseStringPair(const tinyxml2::XMLElement* mom, 
                                             const std::pair<std::string,std::string>& tags);

		std::pair<int,int> ParseIntPair(const tinyxml2::XMLElement* mom, 
                                       const std::pair<std::string,std::string>& tags);

		void ParseParamElement(const tinyxml2::XMLElement* param) const;
		void ParseArrayElement(const tinyxml2::XMLElement* array) const;
		tinyxml2::XMLDocument fDoc;

		std::string fDebugLevel;
		std::string fRootName;
		std::pair<std::string,std::string> fStags;
		std::pair<std::string,std::string> fAlgosPathName;
		std::pair<std::string,std::string> fDstPathName;
		std::pair<std::string,std::string> fHistoPathName;
		std::pair<int,int> fEvents;
		std::string fDebug;

		std::vector<std::string> fAlgoNames;
    std::vector<std::string> fAlgoPath;
			
	};

	typedef SingletonTemplate<AlexConf> AlexConfigure;   // Global declaration

}
#endif