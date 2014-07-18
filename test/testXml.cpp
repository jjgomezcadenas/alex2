// Test for AConfig
// Adapted for NEXT by J.J. Gomez-Cadenas, 2014


#include <alex/AConfigSvc.h>
#include <alex/LogUtil.h>

using namespace tinyxml2;
using namespace alex;

int main(int argc, char **argv)
{
	InitLogger();
  	log4cpp::Category& klog = log4cpp::Category::getRoot();

	alex::AConfigSvc::Instance().ParseParamFile("testXml.xml");
    return 0;
 }