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
using std::stringstream;

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
	AlexConfigure::Instance().Init("DEBUG","AlexConfig");

  if( argc != 1)   
  { 
    cout <<" Usage: xConfigure <pathToAlgos>" << endl;
    return -1;
  }

  stringstream ss;
  ss << argv[0];
  // string pathToMain;
  // ss << argv[1];
  string pathToAlgos;

	string pathToMain ="/Users/jjgomezcadenas/Development/devnext/alex2/ireneMain/";
	//string pathToAlgos ="/Users/jjgomezcadenas/Development/devnext/alex2/irene/";

	string alexConf="AlexConfig.xml";
	string aConfHeader="AConf.hh";
	string aConfCpp="AConf.cxx";
	string algoCpp="AlgoAux.cxx";
	string registerAlgos="RegisterAlgos.cxx";
  string registerAlgosHeader="RegisterAlgosHeader.hh";
	string algoHeader="AlgoHeaders.hh";
  
  klog << log4cpp::Priority::INFO 
        << " path to Algos=" << pathToAlgos;
  klog << log4cpp::Priority::INFO 
        << " path to Main=" << pathToMain;
	
  string pathAlexConf = PathFromStrings(pathToAlgos,alexConf);
  klog << log4cpp::Priority::INFO 
        << " Parse AlexConfig file located at=" << pathAlexConf;

	AlexConfigure::Instance().ParseConfiguration(pathAlexConf);
	
	{
  	string pathAConfHeader = PathFromStrings(pathToAlgos,aConfHeader);
  	klog << log4cpp::Priority::INFO 
    	    << " Write AConf header file at=" << pathAConfHeader;
		std::ofstream out(pathAConfHeader.c_str());
  	out << AlexConfigure::Instance().WriteAConfHeader();
  	out.close();
	}
	{
  	string pathAConfCpp = PathFromStrings(pathToAlgos,aConfCpp);
  	klog << log4cpp::Priority::INFO 
    	    << " Write AConf cpp file at=" << pathAConfCpp;
		std::ofstream out(pathAConfCpp.c_str());
  	out << AlexConfigure::Instance().WriteAConfCPP();
  	out.close();
	}
	{
  	string pathRegisterAlgos = PathFromStrings(pathToAlgos,registerAlgos);
  	klog << log4cpp::Priority::INFO 
    	    << " Write registerAlgo file at=" << pathRegisterAlgos;
		std::ofstream out(pathRegisterAlgos.c_str());
  	out << AlexConfigure::Instance().WriteRegisterAlgos();
  	out.close();
	}
	{
  	string pathAlgoHeader = PathFromStrings(pathToAlgos,algoHeader);
  	klog << log4cpp::Priority::INFO 
    	    << " Write algoHeader file at=" << pathAlgoHeader;
		std::ofstream out(pathAlgoHeader.c_str());
  	out << AlexConfigure::Instance().WriteAlgoHeaderFile();
  	out.close();
	}
  {
    string pathRegisterAlgoHeader = PathFromStrings(pathToAlgos,registerAlgosHeader);
    klog << log4cpp::Priority::INFO 
          << " Write Register algoHeader file at=" << pathRegisterAlgoHeader;
    std::ofstream out(pathRegisterAlgoHeader.c_str());
    out << AlexConfigure::Instance().WriteRegisterAlgosHeader();
    out.close();
  }


	std::vector<std::string> ah = AlexConfigure::Instance().WriteAlgoHeaders();
	size_t i=0;
	for (auto algoName : AlexConfigure::Instance().AlgoNames()) 
	{
		string algoHeader = MergeStrings(algoName,".hh");
  	string pathAlgoHeader = PathFromStrings(pathToAlgos,algoHeader);
  	klog << log4cpp::Priority::INFO 
    	    << " Write Algo header file at=" << pathAlgoHeader;
		std::ofstream out(pathAlgoHeader.c_str());
  	out << ah.at(i);
  	out.close();
  	i++;
	}
	{
  	string pathAlgoCpp = PathFromStrings(pathToAlgos,algoCpp);
  	klog << log4cpp::Priority::INFO 
    	    << " Write Algo cpp file at=" << pathAlgoCpp;
		std::ofstream out(pathAlgoCpp.c_str());
  	out << AlexConfigure::Instance().WriteAlgoCPP() ;
  	out.close();
	}

  return 0;
 }