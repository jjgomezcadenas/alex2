#include "NConstructIBetas.hh"
#include <alex/ISvc.h>
namespace alex {

  //--------------------------------------------------------------------
  bool NConstructIBetas::Init()
  //--------------------------------------------------------------------
  {
    // Init the IreneManager
    ISvc::Instance().Init("DEBUG");
   
    log4cpp::Category& klog = log4cpp::Category::getRoot(); 
    klog << log4cpp::Priority::DEBUG << " fMaxDist = " << fMaxDist;
    return true;
  }
  //--------------------------------------------------------------------
  bool NConstructIBetas::Execute()
  //--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::DEBUG << " NConstructIBetas::Execute";

    // Call the ISvc function to create the IBetas.
    ISvc::Instance().CreateTracks(fMaxDist);
    klog << log4cpp::Priority::DEBUG << " -- Created IBeta objects";

    // Fill a histogram with the number of IBetas for each event.
    std::vector<const IBeta*> iBetas = ISvc::Instance().GetIBetas();
    fH1_nIBetas->Fill(iBetas.size());

    return true;
  }
  //--------------------------------------------------------------------
  bool NConstructIBetas::End()
  //--------------------------------------------------------------------
  {
    return true;
  }
}
