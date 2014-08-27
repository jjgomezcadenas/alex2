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

#include <irene/AlgoHeaders.hh> 
#include <irene/RegisterAlgosHeader.hh>
#include <irene/AConf.hh> // generated

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
  
  //Alex::Instance().InitHistoFile(histoPath);
  TFile* fHistoFile = new TFile(histoPath.c_str(),"RECREATE");
  
	klog << log4cpp::Priority::INFO 
        << " Instantiate and register algos " ;

  // Algos must be initiated before we open the irene TFile, so that their pointers
  // stay in the directory of the histogram file
  //--------

  
  alex::RegisterAlgos();
	alex::Alex::Instance().InitAlgorithms();

  //----
  // Instantiate irene Event
  irene::Event* ievt = new irene::Event();

  //Get path, open DST file, set Tree
  string fp = PathFromStrings(aconf.DstPath(),aconf.DstName());
  klog << log4cpp::Priority::INFO 
        << " Open DST file =" << fp;
  
  TFile* ifile = new TFile(fp.c_str(), "READ");
  TTree* fEvtTree = dynamic_cast<TTree*>(ifile->Get("EVENT"));
  fEvtTree->SetBranchAddress("EventBranch", &ievt);

  klog << log4cpp::Priority::INFO << "number of entries in Irene Tree = " 
  << fEvtTree->GetEntries();
  klog << log4cpp::Priority::INFO << "number of events required = " 
  << aconf.EventsToRun();

  auto nRun = std::min(aconf.EventsToRun(), (int) fEvtTree->GetEntries());
  klog << log4cpp::Priority::INFO << "number of events to run  = " 
  << nRun ;
  
	//-----------Event loop --------
	int nb;
	klog << log4cpp::Priority::INFO 
        << " Start loop " ;
	
  int nev =0;
  int npass = 0;
  int nfail = 0;
  for (int ivt = 0; ivt < nRun; ivt++)
  {
    nb = fEvtTree->GetEntry(ivt);
    //nb = ISvc::Instance().DstGetEntry(ivt);
    ISvc::Instance().LoadEvent(ievt);
    
    if (ivt%aconf.EventsToDebug() ==0)
    {
      klog << log4cpp::Priority::INFO 
          << " read event " << ivt << " nb = " << nb;
    }

    klog << log4cpp::Priority::DEBUG 
        << " Executing algos  " ;

    nev++;
    bool test = alex::Alex::Instance().ExecuteAlgorithms();
    if (test == true)
      npass++;
    else
      nfail++;
    
  }

  klog << log4cpp::Priority::INFO 
        << " Ending...  " ;
  alex::Alex::Instance().EndAlgorithms();
  klog << log4cpp::Priority::INFO  << "Read " 
  << nev << " events" ;
  klog << log4cpp::Priority::INFO  << "Passed selection =" 
  << npass << " Failed selection =" << nfail ;

  // gFile->ls();
  // gFile = fHistoFile;
  // gFile->ls();

  fHistoFile->Write();
  fHistoFile->Close();
  //Alex::Instance().WriteHistoFile();
  //Alex::Instance().CloseHistoFile();
  //Alex::Instance().ClearAlgorithms();

   //------------
    //theApp->Run();
   return 0;
 }