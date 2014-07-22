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

	AlexConfigure::Instance().Init("DEBUG","AlexConfig");

	klog << log4cpp::Priority::INFO 
        << " Parse AlexConfig file " ;

	AlexConfigure::Instance().
	ParseConfiguration("/Users/jjgomezcadenas/Development/devnext/alex2/toy/AlexConfig.xml");

  return 0;
 }