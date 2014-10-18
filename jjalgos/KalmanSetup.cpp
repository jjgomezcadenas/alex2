#include "KalmanSetup.hh"
#include <alex/KFSetup.h>
#include <alex/LogUtil.h>
namespace alex {

//--------------------------------------------------------------------
  bool KalmanSetup::Init()
//--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::DEBUG << " KalmanSetup::Execute" ;

    KFSetup::Instance().Init(gas, model, Pr, B);
    KFSetup::Instance().SetFitParameters (maxChi2, maxOutliers, maxExtrapFailures);
    KFSetup::Instance().SetVerbosity(fitterVerbosity,navigationVerbosity,modelVerbosity,
                                     matchingVerbosity);

    KFSetup::Instance().Print();
  	return true;
  }
//--------------------------------------------------------------------
  bool KalmanSetup::Execute()
//--------------------------------------------------------------------
  {
  	

		return true;
  }
//--------------------------------------------------------------------
  bool KalmanSetup::End()
//--------------------------------------------------------------------
  {
  	return true;
  }
}
