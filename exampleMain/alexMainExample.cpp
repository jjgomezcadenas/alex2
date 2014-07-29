// alexMain example
// J.J. Gomez-Cadenas, August 2014


#include <alex/Alex.h>
#include <alex/StringOperations.h>

#include <TApplication.h>
//#include <Riostream.h>
#include <TVector3.h>
#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <map>

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdio>
#include <cstring>

#include <example/AEx1.hh> // generated
#include <example/AEx2.hh> // generated
#include <example/AConf.hh> // generated
#include <example/ExData.h>


using namespace alex;
using std::string; 
using std::cout; 
using std::endl; 
using std::ostream;
using std::ifstream;
using std::vector;

int main(int argc, char **argv)
{
	// to display histos online
  //TApplication* theApp = new TApplication("App", &argc, argv);

	InitLogger();
  log4cpp::Category& klog = log4cpp::Category::getRoot();
  
  auto aconf = AAConf();
  
  klog << log4cpp::Priority::INFO 
        << " Init Alex with debug level " << aconf.DebugLevel();

  alex::Alex::Instance().Init(aconf.DebugLevel());

  ifstream filein;
  string filePath = PathFromStrings(aconf.DstPath(),aconf.DstName());

  klog << log4cpp::Priority::INFO 
        << " Open data file =" << filePath;

	//filein.open(Form(filePath.c_str()));
  filein.open(filePath);

  string histoPath = PathFromStrings(aconf.HistoPath(),aconf.HistoName());

  klog << log4cpp::Priority::INFO 
        << " Init Histo file =" << histoPath;

	alex::Alex::Instance().InitHistoFile(histoPath);

	klog << log4cpp::Priority::INFO 
        << " Instantiate and register algos " ;

//--------
// This can be auto-generated
	auto aEx1 = new alex::AEx1();
  auto aEx2 = new alex::AEx2();
	
	alex::Alex::Instance().RegisterAlgorithm(aEx1);
  alex::Alex::Instance().RegisterAlgorithm(aEx2);

  //----

	alex::Alex::Instance().InitAlgorithms();

	
	auto nlines =0;
	double x,y,z;

//------------
	klog << log4cpp::Priority::INFO 
        << " Start loop " ;

 
	for(auto i=0; i<= aconf.EventsToRun(); ++i) 
	{
    filein >> x >> y >> z;
    if (!filein.good()) break;
    if (nlines < aconf.EventsToDebug()) printf("x=%8f, y=%8f, z=%8f\n",x,y,z);
    nlines++;
    
    klog << log4cpp::Priority::DEBUG 
        << " instantiating ExData " ;
    auto exData = new ExData();

    klog << log4cpp::Priority::DEBUG 
        << " Setting name " ;

		exData->SetName("ExData");

    klog << log4cpp::Priority::DEBUG 
        << " Setting data " ;

		exData->SetData(TVector3(x,y,z));

    klog << log4cpp::Priority::DEBUG 
        << " Registering data " ;
		alex::Alex::Instance().RegisterData(exData->Name(), exData);

    klog << log4cpp::Priority::DEBUG 
        << " Executing algos  " ;

    alex::Alex::Instance().ExecuteAlgorithms();

    klog << log4cpp::Priority::DEBUG 
        << " Clearing data  " ;
    alex::Alex::Instance().ClearData();
  }

  klog << log4cpp::Priority::INFO 
        << " Ending...  " ;
  alex::Alex::Instance().EndAlgorithms();
  klog << log4cpp::Priority::INFO  << "Read " 
  << nlines << " lines" ;


   alex::Alex::Instance().ClearAlgorithms();
   alex::Alex::Instance().WriteHistoFile();
   alex::Alex::Instance().CloseHistoFile();

   //------------
    //theApp->Run();
   return 0;
 }