// Test for AConfig
// Adapted for NEXT by J.J. Gomez-Cadenas, 2014


#include <alex/Alex.h>

#include <TApplication.h>
#include <Riostream.h>
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
#include <cstdio>
#include <cstring>

#include <alex/AlgoHeaders.h> // generated
#include <alex/AAConf.h> // generated
#include <alex/ExData.h>


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
  TApplication* theApp = new TApplication("App", &argc, argv);

	InitLogger();
  log4cpp::Category& klog = log4cpp::Category::getRoot();

  klog << log4cpp::Priority::INFO 
        << " Init Alex  " ;

	alex::Alex::Instance().Init("DEBUG");

	klog << log4cpp::Priority::INFO 
        << " Open data file " ;

	ifstream filein;
  aconf = AAConf();
  string filePath = PathFromStrings(aconf.DstPath(),aconf.DstName());
	filein.open(Form(filePath.c_str()));

	klog << log4cpp::Priority::INFO 
        << " Init Histo file " ;

  string histoPath = PathFromStrings(aconf.HistoPath(),aconf.HistoName());
	alex::Alex::Instance().InitHistoFile(histoPath.c_str()));

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

	while (1) 
	{
    filein >> x >> y >> z;
    if (!filein.good()) break;
    if (nlines < 5) printf("x=%8f, y=%8f, z=%8f\n",x,y,z);
    nlines++;
    // IData* itoyData = alex::Alex::Instance().RetrieveData("toyData");
    // toyData =dynamic_cast<ToyData*> (itoyData);

    ToyData* toyData = new ToyData();
		toyData->SetName("toyData");
		toyData->SetData(TVector3(x,y,z));
		alex::Alex::Instance().RegisterData(toyData->Name(), toyData);

		if (nlines < 10)
		{
			klog << log4cpp::Priority::DEBUG << " Serialize data " ;
      klog << log4cpp::Priority::DEBUG << toyData->Serialize();

      // klog << log4cpp::Priority::DEBUG 
      // 	  << " Recreate data " ;
      // toyData->Recreate(toyData->Serialize());
    }
    
    alex::Alex::Instance().ExecuteAlgorithms();
    alex::Alex::Instance().ClearData();
  }

   //toyAlgo->End();
  alex::Alex::Instance().EndAlgorithms();
   cout << "Read " << nlines << " lines" << endl;


   alex::Alex::Instance().ClearAlgorithms();
   alex::Alex::Instance().WriteHistoFile();
   alex::Alex::Instance().CloseHistoFile();

   //------------
    theApp->Run();
   return 0;
 }