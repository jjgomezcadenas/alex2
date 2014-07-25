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

#include <alex/ToyAnalysis.h>
#include <alex/ToyData.h>

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
	filein.open(Form("basic.dat"));

	klog << log4cpp::Priority::INFO 
        << " Init Histo file " ;

	alex::Alex::Instance().InitHistoFile("toyHistos.root");

	klog << log4cpp::Priority::INFO 
        << " Instantiate and register algos " ;

	auto toyAlgo = new alex::ToyAnalysis();
	toyAlgo->SetName("toyAlgo");
	alex::Alex::Instance().RegisterAlgorithm(toyAlgo);

	alex::Alex::Instance().InitAlgorithms();

	// ToyData* toyData = new ToyData();
	// toyData->SetName("toyData");
	// alex::Alex::Instance().RegisterData(toyData->Name(), toyData);

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