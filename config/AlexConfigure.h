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

	private:
		tinyxml2::XMLElement* ParseRoot() ;
		std::pair<std::string,std::string> ParseStringPair(const tinyxml2::XMLElement* mom, 
                                             const std::pair<std::string,std::string>& tags);

		std::pair<int,int> ParseIntPair(const tinyxml2::XMLElement* mom, 
                                       const std::pair<std::string,std::string>& tags);
		tinyxml2::XMLDocument fDoc;

		std::string fDebugLevel;
		std::string fRootName;
			
	};

	typedef SingletonTemplate<AlexConf> AlexConfigure;   // Global declaration

}
#endif