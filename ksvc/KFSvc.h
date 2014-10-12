/* -*- mode: c++ -*- */
#ifndef KFSVC_
#define KFSVC_
/*
 KFSvc: Manages Kalman Filter
 

 JJGC, October, 2014.
*/
#include <alex/SingletonTemplate.h>

#include <alex/KFRecpack.h>

#include <string>
#include <vector>
#include <utility>
#include <memory>



namespace alex {

  class Hit;
    
  class KFSvcManager {
	
    public:
		  KFSvcManager(){};
      virtual ~KFSvcManager(){};
      void Init();
		  
      RP::Trajectory CreateTrajectory(std::vector<const Hit* > hits, 
                                      std::vector<double> hitErrors) ;

      RP::State SeedState(std::vector<double> v0, std::vector<double> p0) ;

      bool FitTrajectory(RP::Trajectory traj,RP::State seed);
     
      std::string Model() const;
      int ModelDim() const ;

      double X0() const {return fX0;}
      double dEdX() const {return fDedx;}
      RP::EVector BField() const {return fBField;}
      RP::EMatrix MeasurementCovariance() const {return fCov; }
      RP::measurement_vector MeasurementVector() const {return fMeas;}
      
 			
  //   protected:

      void InitializeManagerGeometry();
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
			
	};
  typedef alex::SingletonTemplate<KFSvcManager> KFSvc; 
}
#endif