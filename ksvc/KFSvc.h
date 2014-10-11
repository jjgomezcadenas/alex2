#ifndef KFSVC_
#define KFSVC_
/*
 KFSvc: Manages Kalman Filter
 

 JJGC, April, 2014.
*/

#include <recpack/RecPackManager.h>
#include <recpack/string_tools.h>
#include <recpack/stc_tools.h>
#include <recpack/ERandom.h>
#include <recpack/HelixEquation.h>
#include <recpack/LogSpiralEquation.h>
#include <recpack/ParticleState.h>
#include <recpack/RayTool.h>
#include <recpack/EGeo.h>
#include <recpack/Definitions.h>
#include <recpack/Trajectory.h>
#include <recpack/KalmanFitter.h>
#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <alex/SingletonTemplate.h>



namespace alex {
    
  class KFSvcManager {
	
    public:
		  KFSvcManager(){};
      virtual ~KFSvcManager(){};
      void Init();
		  bool Fit();
      RP::Trajectory CreateTrajectory() ;
      void SeedState(const Trajectory& traj, State& state) ;
    
      std::string Model() const;
      int ModelDim() const ;
      
      
 			
    protected:

      void CreateSetup();
      void SetVerbosity();

      std::string fModel;        // model for this setup 
      int fDim;                  // dimension for this setup 

      double fX0;                // radiation length for this setup
      double fDedx;              // dedx for this setup
    
		  RP::RecPackManager fRPMan;   //Recpack magnager
      RP::Setup fSetup;                // RP setup

  // volume and surface properties cannot be local variables 
  // because RecPack takes a reference to the property

      RP::EVector fBField;         // (Bx,By,Bz)
      RP::EMatrix fCov;
      RP::EVector xaxis;
      RP::EVector yaxis;
      RP::EVector zaxis;
  
      RP::measurement_vector fMeas;
      RP::Trajectory fSimTraj;
			
	};
  typedef SingletonTemplate<KFSetupManager> KFSvc; 
}
#endif