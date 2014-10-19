#include "NConstructRBetas.hh"
#include <alex/ISvc.h>
namespace alex {

  //--------------------------------------------------------------------
  bool NConstructRBetas::Init()
  //--------------------------------------------------------------------
  {
    // Init the IreneManager
    ISvc::Instance().Init("DEBUG");
   
    log4cpp::Category& klog = log4cpp::Category::getRoot(); 
    klog << log4cpp::Priority::DEBUG << " fCoreDist = " << fCoreDist;
    klog << log4cpp::Priority::DEBUG << " fSparseWidth = " << fSparseWidth;
    klog << log4cpp::Priority::DEBUG << " fBlobRadius = " << fBlobRadius;

    return true;
  }
  //--------------------------------------------------------------------
  bool NConstructRBetas::Execute()
  //--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::DEBUG << " NConstructRBetas::Execute";

    // Call the ISvc function to create the RBeta.
    ISvc::Instance().CreateRTracks(fCoreDist,fSparseWidth,fBlobRadius);
    klog << log4cpp::Priority::DEBUG << " -- Created RBeta object";

    return true;
  }
  //--------------------------------------------------------------------
  bool NConstructRBetas::End()
  //--------------------------------------------------------------------
  {
    return true;
  }
}
