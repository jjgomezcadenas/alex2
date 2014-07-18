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

	private:
		std::vector<IAlgorithm*> fIAlgo;
			
	};

	typedef SingletonTemplate<AlexManager> Alex;   // Global declaration

}
#endif