// Test for AConfig
// Adapted for NEXT by J.J. Gomez-Cadenas, 2014


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

#include <alex/AlexConfigure.h>
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

  klog << log4cpp::Priority::INFO 
        << " Init textAlexConfigure  " ;

	AlexConfigure::Instance().Init("INFO","AlexConfig");

	klog << log4cpp::Priority::INFO 
        << " Parse AlexConfig file " ;

	AlexConfigure::Instance().
	ParseConfiguration("/Users/jjgomezcadenas/Development/devnext/alex2/example/AlexConfig.xml");
	
	
	klog << log4cpp::Priority::DEBUG << "SerializeAlgoNames";
	klog << log4cpp::Priority::DEBUG << AlexConfigure::Instance().SerializeAlgoNames();
	klog << log4cpp::Priority::DEBUG << "SerializeAlgoPaths";
	klog << log4cpp::Priority::DEBUG << AlexConfigure::Instance().SerializeAlgoPaths();
	klog << log4cpp::Priority::DEBUG << "SerializeAlgoParam";
	klog << log4cpp::Priority::DEBUG << AlexConfigure::Instance().SerializeAlgoParam();
	klog << log4cpp::Priority::DEBUG << "SerializeAlgoArray";
	klog << log4cpp::Priority::DEBUG << AlexConfigure::Instance().SerializeAlgoArray();
	klog << log4cpp::Priority::DEBUG << "SerializeAlgoH1D";
	klog << log4cpp::Priority::DEBUG << AlexConfigure::Instance().SerializeAlgoH1D();
	klog << log4cpp::Priority::DEBUG << "SerializeAlgoH2D";
	klog << log4cpp::Priority::DEBUG << AlexConfigure::Instance().SerializeAlgoH2D();

	klog << log4cpp::Priority::INFO << "Write AConf Header";
	klog << log4cpp::Priority::INFO << AlexConfigure::Instance().WriteAConfHeader();
	klog << log4cpp::Priority::INFO << "Write AConf CPP";
	klog << log4cpp::Priority::INFO << AlexConfigure::Instance().WriteAConfCPP();
	klog << log4cpp::Priority::INFO << "Write Algo Header";
	klog << log4cpp::Priority::INFO << AlexConfigure::Instance().WriteAlgoHeader();
	klog << log4cpp::Priority::INFO << "Write Algo CPP";
	klog << log4cpp::Priority::INFO << AlexConfigure::Instance().WriteAlgoCPP();
	
  return 0;
 }