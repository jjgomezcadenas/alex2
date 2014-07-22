// Test for AConfig
// Adapted for NEXT by J.J. Gomez-Cadenas, 2014


#include <alex/Alex.h>

using namespace alex;

int main(int argc, char **argv)
{
	InitLogger();
  	log4cpp::Category& klog = log4cpp::Category::getRoot();

  	klog << log4cpp::Priority::INFO 
        << " Init Alex  " ;

	alex::Alex::Instance().Init("DEBUG");
    return 0;
 }