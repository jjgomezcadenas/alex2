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

class AXml {
	public:
		AXml(){};
		virtual ~AXml(){};
		void ParseParamFile(std::string paramFile);

	private:
		void ParseParamElement(const tinyxml2::XMLElement* param) const;
		void ParseArrayElement(const tinyxml2::XMLElement* array) const;
		tinyxml2::XMLDocument fDoc;
			
	};

	typedef SingletonTemplate<AXml> AConfigSvc;   // Global declaration

}
#endif