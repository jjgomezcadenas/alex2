// ireneMain example
// J.J. Gomez-Cadenas, August 2014


#include <alex/Alex.h>
#include <alex/ISvc.h>
#include <alex/StringOperations.h>

#include <TApplication.h>
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

#include <irene/IElectrons.hh> // generated
#include <irene/AConf.hh> // generated
// #include <example/AEx2.hh> // generated

// #include <example/ExData.h>


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
  //histogram file;
  string histoPath = PathFromStrings(aconf.HistoPath(),aconf.HistoName());
  klog << log4cpp::Priority::INFO 
        << " Init Histo file =" << histoPath;
  Alex::Instance().InitHistoFile(histoPath);

  // Instantiate irene Event
  irene::Event* ievt = new irene::Event();

  //Get path
  string fp = PathFromStrings(aconf.DstPath(),aconf.DstName());
  klog << log4cpp::Priority::INFO 
        << " Open data file =" << fp;

  //ISvc::Instance().InitDst(fp,ievt);
  TFile* ifile = new TFile(fp.c_str(), "READ");
  TTree* fEvtTree = dynamic_cast<TTree*>(ifile->Get("EVENT"));
  fEvtTree->SetBranchAddress("EventBranch", &ievt);

  // klog << log4cpp::Priority::INFO << "number of entries in Irene Tree = " 
  // << ISvc::Instance().DstEntries();
  klog << log4cpp::Priority::INFO << "number of entries in Irene Tree = " 
  << fEvtTree->GetEntries();
  klog << log4cpp::Priority::INFO << "number of events required = " 
  << aconf.EventsToRun();

  //auto nRun = std::min(aconf.EventsToRun(), ISvc::Instance().DstEntries());
  auto nRun = std::min(aconf.EventsToRun(), (int) fEvtTree->GetEntries());
  klog << log4cpp::Priority::INFO << "number of events to run  = " 
  << nRun ;
  
	klog << log4cpp::Priority::INFO 
        << " Instantiate and register algos " ;

//--------
// This can be auto-generated
	auto iElectrons = new alex::IElectrons();
  //auto aEx2 = new alex::AEx2();
	
	Alex::Instance().RegisterAlgorithm(iElectrons);
  //alex::Alex::Instance().RegisterAlgorithm(aEx2);

  //----

	alex::Alex::Instance().InitAlgorithms();

	//-----------Event loop --------
	int nb;
	klog << log4cpp::Priority::INFO 
        << " Start loop " ;
	
  int nev =0;
  for (int ivt = 0; ivt < nRun; ivt++)
  {
    nb = fEvtTree->GetEntry();
    //nb = ISvc::Instance().DstGetEntry(ivt);
    ISvc::Instance().LoadEvent(ievt);
    
    klog << log4cpp::Priority::DEBUG 
        << " Executing algos  " ;
    alex::Alex::Instance().ExecuteAlgorithms();
    nev++;
  }

  klog << log4cpp::Priority::INFO 
        << " Ending...  " ;
  alex::Alex::Instance().EndAlgorithms();
  klog << log4cpp::Priority::INFO  << "Read " 
  << nev << " events" ;

   alex::Alex::Instance().ClearAlgorithms();
   alex::Alex::Instance().WriteHistoFile();
   alex::Alex::Instance().CloseHistoFile();

   //------------
    //theApp->Run();
   return 0;
 }