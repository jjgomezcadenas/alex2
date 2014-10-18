#include "KalmanSetup.hh"
#include <alex/KFSetup.h>
#include <alex/KFSvc.h>
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

    KFSvc::Instance().Init();  // this must be called after KFSetup
  	return true;
  }
//--------------------------------------------------------------------
  bool KalmanSetup::Execute()
//--------------------------------------------------------------------
  {
  	
    // you must call here the methods to get a trajectory and fit it:
    /*
    Trajectory* trj = CreateTrajectory(std::vector<const Hit* > hits, 
                                      std::vector<double> hitErrors) ;

    State* seed =SeedState(std::vector<double> v0, std::vector<double> p0) ;

    bool status = FitTrajectory(*traj, *seed);

    // Here we must access the Trajectory information and retrieve stuff. 

    This is Sasha's e-mail, don't have much more info. 

    Regarding which info can  be retrieved from nodes as a result of the
kalman filtering:
different state hypervector values and residuals are stored with the
corresponding names:
E.g. to retrieve predicted dynamic vector and matrix from prediction
step:
// get the state
State& state = node.state();
const HyperVector& hv_P = state.hv(RP::predicted);
const EVector& x_P = hv_P.vector();
const EMatrix& C_P = hv_P.matrix();
// Retrieve predicted residual
const HyperVector resHV_P = node.residuals().hv(predicted_name); const
EVector& r_P = resHV_P.vector();
const EMatrix& R_P = resHV_P.matrix();
node.quality(predicted_name,chi2); //or resHV_P.chi2()
available names:
RP:predicted, RP::filtered
RP::smoothed (default)


    delete trj;
    delete seed;
    */

		return true;
  }
//--------------------------------------------------------------------
  bool KalmanSetup::End()
//--------------------------------------------------------------------
  {
  	return true;
  }
}
