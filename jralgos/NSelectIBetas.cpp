#include "NSelectIBetas.hh"
#include <alex/ISvc.h>
namespace alex {

  //--------------------------------------------------------------------
  bool NSelectIBetas::Init()
  //--------------------------------------------------------------------
  {
    // Init the IreneManager
    ISvc::Instance().Init("DEBUG");

    log4cpp::Category& klog = log4cpp::Category::getRoot();    
    klog << log4cpp::Priority::DEBUG << " fSelect = " << fSelect;
    klog << log4cpp::Priority::DEBUG << " fMaxIBetas = " << fMaxIBetas;

    return true;
  }
  //--------------------------------------------------------------------
  bool NSelectIBetas::Execute()
  //--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::DEBUG << " NSelectIBetas::Execute";

    // Perform the IBetas cut if enabled.
    if(fSelect == 1) {
      
      // Get the number of IBeta objects in the event and compare to the max allowed.
      std::vector<const IBeta*> iBetas = ISvc::Instance().GetIBetas();
      if(fMaxIBetas < (int) iBetas.size()) {
        klog << log4cpp::Priority::DEBUG << " -- Event failed IBeta selection cut with " 
             << iBetas.size() << " IBeta objects.";
        return false;
      }
    }

    return true;
  }
  //--------------------------------------------------------------------
  bool NSelectIBetas::End()
  //--------------------------------------------------------------------
  {
    return true;
  }
}
