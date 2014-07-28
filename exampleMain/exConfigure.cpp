// exConfigure
// Configures the setup to run the algos in folder example
// August, 2014


#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <map>

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstring>
#include <fstream>

#include <alex/AlexConfigure.h>
#include <alex/StringOperations.h>
#include <alex/LogUtil.h>

using namespace alex;
using std::string; 
using std::cout; 
using std::endl; 
using std::ostream;
using std::ifstream;
using std::vector;

int main(int argc, char **argv)
{
	
	InitLogger();
	log4cpp::Category& klog = log4cpp::Category::getRoot();
	AlexConfigure::Instance().Init("INFO","AlexConfig");

	string pathToAux ="/Users/jjgomezcadenas/Development/devnext/alex2/example/";
	string pathToAlgos ="/Users/jjgomezcadenas/Development/devnext/alex2/example/";
	string alexConf="AlexConfig.xml";
	string aConfHeader="AConf.hh";
	string aConfCpp="AConf.cxx";
	string algoHeader="AExample.h";
	string algoCpp="AExample.cxx";
  
  klog << log4cpp::Priority::INFO 
        << " path to Algos  " << pathToAlgos;
  klog << log4cpp::Priority::INFO 
        << " path to Aux  " << pathToAux;
	
  string pathAlexConf = PathFromStrings(pathToAlgos,alexConf);
  klog << log4cpp::Priority::INFO 
        << " Parse AlexConfig file located at " << pathAlexConf;

	AlexConfigure::Instance().ParseConfiguration(pathAlexConf);
	
	klog << log4cpp::Priority::INFO << "Write AConf Header";
	klog << log4cpp::Priority::INFO << AlexConfigure::Instance().WriteAConfHeader();
	klog << log4cpp::Priority::INFO << "Write AConf CPP";
	klog << log4cpp::Priority::INFO << AlexConfigure::Instance().WriteAConfCPP();
	klog << log4cpp::Priority::INFO << "Write Algo Header";
	klog << log4cpp::Priority::INFO << AlexConfigure::Instance().WriteAlgoHeader();
	klog << log4cpp::Priority::INFO << "Write Algo CPP";
	klog << log4cpp::Priority::INFO << AlexConfigure::Instance().WriteAlgoCPP();

	{
  	string pathAConfHeader = PathFromStrings(pathToAux,aConfHeader);
  	klog << log4cpp::Priority::INFO 
    	    << " Write AConf header file at " << pathAConfHeader;
		std::ofstream out(pathAConfHeader.c_str());
  	out << AlexConfigure::Instance().WriteAConfHeader();
  	out.close();
	}
	{
  	string pathAConfCpp = PathFromStrings(pathToAux,aConfCpp);
  	klog << log4cpp::Priority::INFO 
    	    << " Write AConf cpp file at " << pathAConfCpp;
		std::ofstream out(pathAConfCpp.c_str());
  	out << AlexConfigure::Instance().WriteAConfCPP();
  	out.close();
	}

	std::vector<std::string> ah = AlexConfigure::Instance().WriteAlgoHeaders();
	size_t i=0;
	for (auto algoName : AlexConfigure::Instance().AlgoNames()) 
	{
		string algoHeader = MergeStrings(algoName,".hh");
  	string pathAlgoHeader = PathFromStrings(pathToAux,algoHeader);
  	klog << log4cpp::Priority::INFO 
    	    << " Write Algo header file at " << pathAlgoHeader;
		std::ofstream out(pathAlgoHeader.c_str());
  	out << ah.at(i);
  	out.close();
  	i++;
	}
	{
  	string pathAlgoCpp = PathFromStrings(pathToAux,algoCpp);
  	klog << log4cpp::Priority::INFO 
    	    << " Write Algo cpp file at " << pathAlgoCpp;
		std::ofstream out(pathAlgoCpp.c_str());
  	out << AlexConfigure::Instance().WriteAlgoCPP() ;
  	out.close();
	}

  return 0;
 }