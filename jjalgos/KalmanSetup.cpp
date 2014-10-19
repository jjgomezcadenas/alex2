#include "KalmanSetup.hh"
#include <alex/KFSetup.h>
#include <alex/KFSvc.h>
#include <alex/ISvc.h>
#include <alex/LogUtil.h>

using namespace std;
namespace alex {

//--------------------------------------------------------------------
  bool KalmanSetup::Init()
//--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::DEBUG << " KalmanSetup::Init" ;

    KFSetup::Instance().Init(gas, model, Pr, B);
    KFSetup::Instance().SetFitParameters (maxChi2, maxOutliers, maxExtrapFailures);
    KFSetup::Instance().SetVerbosity(fitterVerbosity,navigationVerbosity,modelVerbosity,
                                     matchingVerbosity);


    klog << log4cpp::Priority::DEBUG << " KFSetup initialized" ;
    KFSetup::Instance().Print();

    string input ;
    cout << "Return to continue:\n>";
    getline(cin, input);

    KFSvc::Instance().Init();  // this must be called after KFSetup

    klog << log4cpp::Priority::DEBUG << " KFSvc initialized" ;
    KFSetup::Instance().Print();

  	return true;
  }
//--------------------------------------------------------------------
  bool KalmanSetup::Execute()
//--------------------------------------------------------------------
  {
  	
   log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::DEBUG << " KalmanSetup::Execute" ;

    RBeta* rBeta = ISvc::Instance().GetRBeta();
    std::vector<const Hit*> effHits = rBeta->GetEffHits();

    ISvc::Instance().CreateKFObjects(QMAX,1,true,true);

    std::vector<double> hitErrors = ISvc::Instance().GetKFMErrors();
    std::vector<double> v0 = ISvc::Instance().GetKFv0();
    std::vector<double> p0 = ISvc::Instance().GetKFp0();

    string input ;

    klog << log4cpp::Priority::DEBUG << " Creating trajectory" ;
    getline(cin, input);
   
    Trajectory* trj = KFSvc::Instance().CreateTrajectory(effHits, hitErrors);
                                      
    klog << log4cpp::Priority::DEBUG << " Creating seed" ;
    getline(cin, input);

    State* seed =KFSvc::Instance().SeedState(v0, p0) ;

   
    klog << log4cpp::Priority::DEBUG << " Fitting trajectory" ;
    getline(cin, input);

    bool status = KFSvc::Instance().FitTrajectory(*trj, *seed);

    klog << log4cpp::Priority::DEBUG << " Fitting trajectory, result " << status;
    getline(cin, input);

    klog << log4cpp::Priority::DEBUG << " trajectory length = " << trj->length();

    // Get the list of nodes:

    klog << log4cpp::Priority::DEBUG << " Get the list of nodes and loop over them " ;

    const std::vector<Node*> tnodes =  trj->nodes();

    auto inode =0;
    for (auto node: tnodes)
    {
      // get the state X and C

      State state = node->state();
      const HyperVector hv_P = state.hv(RP::predicted);  //predicted
      const EVector x_P = hv_P.vector();  // vector state
      const EMatrix C_P = hv_P.matrix();  //cov matrix

      // Retrieve predicted residual
      /*
      available names:
      RP:predicted, RP::filtered
      RP::smoothed (default)
      */
      HyperVector resHV_P = node->residuals().hv(RP::predicted); 
      const EVector r_P = resHV_P.vector();
      const EMatrix R_P = resHV_P.matrix();
      double tchi2 = resHV_P.chi2();

      klog << log4cpp::Priority::DEBUG << " For node " << inode;
      klog << log4cpp::Priority::DEBUG 
      << " X_P[0] " << x_P[0]
      << " X_P[1] " << x_P[1]
      << " X_P[2] " << x_P[2]
      << " X_P[3] " << x_P[3]
      << " X_P[4] " << x_P[4];

      klog << log4cpp::Priority::DEBUG<< " chi2 = " << tchi2 ;

      inode++;
    }

    getline(cin, input);
    delete trj;
    delete seed;

		return true;
  }
//--------------------------------------------------------------------
  bool KalmanSetup::End()
//--------------------------------------------------------------------
  {
  	return true;
  }
}
