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

	alex::Alex::Instance().Init("DEBUG");

	ifstream filein;
	filein.open(Form("basic.dat"));

	auto toyAlgo = new alex::ToyAnalysis();
	toyAlgo->SetName("toyAlgo");

	alex::Alex::Instance().RegisterAlgorithm(toyAlgo);

	alex::Alex::Instance().InitAlgorithms();

	ToyData* toyData = new ToyData();
	toyData->SetName("toyData");
	alex::Alex::Instance().RegisterData(toyData->Name(), toyData);

	auto nlines =0;
	double x,y,z;

//------------
	while (1) 
	{
    filein >> x >> y >> z;
    if (!filein.good()) break;
    if (nlines < 5) printf("x=%8f, y=%8f, z=%8f\n",x,y,z);
    nlines++;
    IData* itoyData = alex::Alex::Instance().RetrieveData("toyData");
    toyData =dynamic_cast<ToyData*> (itoyData);
    toyData->SetData(TVector3(x,y,z));
    alex::Alex::Instance().ExecuteAlgorithms();
    //toyAlgo->Execute();
  }

   //toyAlgo->End();
  alex::Alex::Instance().EndAlgorithms();
   cout << "Read " << nlines << " lines" << endl;


  //------------
    theApp->Run();

   alex::Alex::Instance().ClearAlgorithms();
   return 0;
 }