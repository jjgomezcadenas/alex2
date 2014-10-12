/* -*- mode: c++ -*- */
#ifndef _kfsetup_rp___
#define _kfsetup_rp___

#include <string>
#include <sstream>
#include <vector>
#include <CLHEP/Units/SystemOfUnits.h>
#include <alex/SingletonTemplate.h>

using std::string; 
using std::vector;
using std::ostream;

namespace alex {
  class KFSetupManager{
    public:
      KFSetupManager(){};
      virtual ~KFSetupManager(){};
      void Init(string gas, string model, double Pr, double B);
      void SetFitParameters (double maxChi2,double maxOutliers, double maxExtrapFailures,
                            double minDistanceNodeOrdering, double minGoodNodeOrdering  );
      void SetVerbosity(string fitterVerbosity,string navigationVerbosity,string modelVerbosity,
                 string matchingVerbosity);

    private:
      string fGas;  // gas type (Xe, H2, etc)
      string fModel; // helix, sline, etc.

      //verbosity levels
      string fFitterVerbosity;
      string fNavigationVerbosity;
      string fModelVerbosity;
      string fMatchingVerbosity;

      double fPr;         // gas pressure
      double fB;         //magnetic field
      double fMaxChi2 ;  // max local chi2/ndf allowed in fit
      double fMaxOutliers ;  // max number of outliers
      double fMaxExtrapFailures; //max number of extrap failures
      double fMinDistanceNodeOrdering ;// minimum distance between nodes to check the node ordering. 
      double fMinGoodNodeOrdering ; // minimum number of nodes to be check for good ordering
      double fRho1B ; //density of gas at 1 bar
      double fDeDx1B ;      // ionization (mip) at 1 bar
      double fX01B ;    // radiation lenght at 1 bar

      double fRho ; // densitiy at gas pressure
      double fX0 ;        // rad length at gas pressure
      double fDeDx;  //dedx  at gas pressure
    
    public: 
    
      string Gas()const {return fGas;}
      string Model()const {return fModel;}
      double Pressure() const {return fPr;}
      double B() const {return fB;}
    
      double Density() const {return fRho;}
      double X0() const {return fX0;}
      double DEDX() const {return fDeDx;}
      vector<string> FiNaMoMaVerbosity() const;  //FitterNavigationModelMatching

      double MaxChi2()const {return fMaxChi2;}
      double MaxOutliers()const {return fMaxOutliers;}
      double MaxExtrapFailures()const {return fMaxExtrapFailures;}
      double MinDistanceNodeOrdering()const {return fMinDistanceNodeOrdering;}
      double MinGoodNodeOrdering()const {return fMinGoodNodeOrdering;}

      void Info(ostream& s) const;
      void Print() const;
  };
  
  ostream& operator << (ostream& s, const KFSetupManager& p);
  typedef alex::SingletonTemplate<KFSetupManager> KFSetup;   // Global declaration
}

#endif
