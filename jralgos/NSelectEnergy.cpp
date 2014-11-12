#include "NSelectEnergy.hh"
#include <alex/ISvc.h>
namespace alex {

  //--------------------------------------------------------------------
  bool NSelectEnergy::Init()
  //--------------------------------------------------------------------
  {
    // Init the IreneManager
    ISvc::Instance().Init("DEBUG");

    log4cpp::Category& klog = log4cpp::Category::getRoot();    
    klog << log4cpp::Priority::DEBUG << " fSelect = " << fSelect;
    klog << log4cpp::Priority::DEBUG << " fMinEnergy = " << fMinEnergy;

    return true;
  }
  //--------------------------------------------------------------------
  bool NSelectEnergy::Execute()
  //--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::DEBUG << " NSelectEnergy::Execute";

    // Perform the energy cut if enabled.
    if(fSelect == 1) {
      
      // Get the energy in the event and compare to the minimum allowed.
      double etot = 0.;
      IHits truehits = ISvc::Instance().GetTrueHits();
      for(int h = 0; h < (int) truehits.size(); h++) {
        etot += truehits[h].second;
      }
      if(etot < fMinEnergy) {
        klog << log4cpp::Priority::DEBUG << " -- Event failed energy selection cut with energy " 
             << etot << " MeV.";
        return false;
      }
    }

    return true;
  }
  //--------------------------------------------------------------------
  bool NSelectEnergy::End()
  //--------------------------------------------------------------------
  {
    return true;
  }
}
