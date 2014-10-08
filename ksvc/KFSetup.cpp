// A simple class to define a Kalman Filter Setup
// JJ, April 2014


#include <alex/KFSetup.h>
#include <alex/LogUtil.h>
#include <cstdlib>

using namespace CLHEP;

namespace alex {

  void KFSetupManager::Init(string gas, string model, double Pr, double B)
    {
      log4cpp::Category& klog = log4cpp::Category::getRoot();
      fGas = gas;
      fModel = model ;
      fPr = Pr;
      fB = B*tesla;
    
      if (fGas == "Xe")
      {
        fRho1B = 5.48e-3*gram/cm3; //density of xenon at 1 bar
        fDeDx1B = 6.88*keV/cm;      // ionization (mip) at 1 bar
        fX01B = 1547*cm;    // radiation lenght at 1 bar
      }
      else if (fGas == "H2")
      {
        fRho1B = 8.38E-5*gram/cm3; //density of xenon at 1 bar
        fDeDx1B = 3.44E-4*MeV/cm;      // ionization (mip) at 1 bar
        fX01B = 7.5e+5*cm;    // radiation lenght at 1 bar
      }
      else
      {
        std::cout << "ERROR: gas unknown" << std::endl;
            exit (EXIT_FAILURE);
      }

      string s="Gas type = %s Gas density: rho (1bar) =%7.5f";
      s+="gram/cm3 dE/dx =%7.2f keV/cm  X0= %7.1f cm \n";
      klog.info(s.c_str(),
                fGas.c_str(),fRho1B/(gram/cm3),fDeDx1B/(keV/cm), fX01B/cm);

      fRho = fPr*fRho1B; // in good approx. density is proportional to pressure.
      fX0 = fX01B/fPr;        // rad length is prop to inverse of pressure
      fDeDx = fPr*fDeDx1B;  //dedX is proportional to pressure

      s="Gas pressure: =%7.1f bar Gas density: (at P) =%7.5f";
      s+="gram/cm3 dE/dx =%7.1f keV/cm  X0= %7.1f cm \n";
      klog.info(s.c_str(),
                fPr, fRho/(gram/cm3),fDeDx/(keV/cm), fX0/cm);

      klog.info("Model =%s, B (tesla) = %7.1f \n", fModel.c_str(),fB/tesla);
   }

  void KFSetupManager::SetFitParameters (double maxChi2,double maxOutliers, double maxExtrapFailures,
                                  double minDistanceNodeOrdering, double minGoodNodeOrdering)
  {
    fMaxChi2 = maxChi2;
    fMaxOutliers = maxOutliers;  
    fMaxExtrapFailures = maxExtrapFailures; 
    fMinDistanceNodeOrdering = minDistanceNodeOrdering; 
    fMinGoodNodeOrdering = minGoodNodeOrdering;

  }
  void KFSetupManager::SetVerbosity(string fitterVerbosity,string navigationVerbosity,string modelVerbosity,
                 string matchingVerbosity)
  {
    fFitterVerbosity = fitterVerbosity;
    fNavigationVerbosity = navigationVerbosity;
    fModelVerbosity = modelVerbosity;
    fMatchingVerbosity = matchingVerbosity;
  }

  vector<string> KFSetupManager::FiNaMoMaVerbosity() const
  {
    vector<string> vb;
    vb.push_back(fFitterVerbosity);
    vb.push_back(fNavigationVerbosity);
    vb.push_back(fModelVerbosity);
    vb.push_back(fMatchingVerbosity);
    return vb;
  }
  void KFSetupManager::Info(ostream& s) const
  {
    s << std::endl;
    s<<"{KFSetupManager: ";
    s << "Gas=" << fGas << " Pressure(bar) =" 
      << fPr  << " B (tesla) =" << fB/tesla << std::endl;
    s << "X0 (cm) =" << fX0/cm << " dE/dx (keV/cm) = "
     << fDeDx/(keV/cm) << std::endl;
    s << "Max chi2 = " << fMaxChi2 << "}";
    s << std::endl;
    s<<"{KFSetupManager Verbosity: ";
    s << "fitterVerbosity = " << fFitterVerbosity << std::endl;
    s << "navigationVerbosity = " << fNavigationVerbosity << std::endl;
    s << "modelVerbosity = " << fModelVerbosity << std::endl;
    s << "matchingVerbosity = " << fMatchingVerbosity << std::endl;
    s << "}";

  }

  ostream& operator << (ostream& s, const KFSetupManager& p) 
  {
    p.Info(s);
    return s; 
  }
}
