// Generated by AlexConf: do not edit
#include "AHits.hh"
#include <alex/ISvc.h>
namespace alex {

//--------------------------------------------------------------------
  bool AHits::Init()
//--------------------------------------------------------------------
  {
  	return true;
  }
//--------------------------------------------------------------------
  bool AHits::Execute()
//--------------------------------------------------------------------
  {
  	log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::DEBUG << " AHits::Execute" ;

  	IHits vHits = ISvc::Instance().GetTrueHits();
    klog << log4cpp::Priority::DEBUG << " true hits = " << vHits.size();
    
  	fH1_nth->Fill(vHits.size());

		return true;
  }
//--------------------------------------------------------------------
  bool AHits::End()
//--------------------------------------------------------------------
  {
  	return true;
  }
}