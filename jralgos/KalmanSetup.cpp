#include "KalmanSetup.hh"
#include <fstream>
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
    //getline(cin, input);

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
    //std::vector<const Hit*> effHits = rBeta->GetEffHits();

    bool forwardFit = true;
    if(fForwardFit == 0) forwardFit = false;

    ISvc::Instance().CreateKFObjects(QMAX+MELEC,1.,true,forwardFit);

    std::vector<const Hit*> effHits = ISvc::Instance().GetKFHits();
    std::vector<double> hitErrors = ISvc::Instance().GetKFMErrors();
    std::vector<double> v0 = ISvc::Instance().GetKFv0();
    std::vector<double> p0 = ISvc::Instance().GetKFp0();

    // Add the obtained effective hits to a histogram.
    for(int i = 0; i < (int) effHits.size(); i++) {
      double x = effHits[i]->XYZ().X();
      double y = effHits[i]->XYZ().Y();
      fH2_xyMeasured->Fill(x,y,i);
    }

    string input ;

    klog << log4cpp::Priority::DEBUG << " Creating trajectory" ;
    //getline(cin, input);
   
    Trajectory* trj = KFSvc::Instance().CreateTrajectory(effHits, hitErrors);
                                      
    klog << log4cpp::Priority::DEBUG << " Creating seed" ;
    //getline(cin, input);

    State* seed =KFSvc::Instance().SeedState(v0, p0) ;

   
    klog << log4cpp::Priority::DEBUG << " Fitting trajectory" ;
    //getline(cin, input);

    bool status = KFSvc::Instance().FitTrajectory(*trj, *seed);

    klog << log4cpp::Priority::DEBUG << " Fitting trajectory, result " << status;
    //getline(cin, input);

    klog << log4cpp::Priority::DEBUG << " trajectory length = " << trj->length();

    // Get the list of nodes:

    klog << log4cpp::Priority::DEBUG << " Get the list of nodes and loop over them " ;

    // Print a file analogous to that for a ToyMC particle.
    char tempStr[200];

    if(forwardFit) {
      sprintf(tempStr,"mkdir -p %s/%s/toyMC/%s",outPath.c_str(),runName.c_str(),runName.c_str());
      system(tempStr);
      sprintf(tempStr,"%s/%s/toyMC/%s/%s_%i.dat",outPath.c_str(),runName.c_str(),runName.c_str(),runName.c_str(),fEvent);
    }
    else {
      sprintf(tempStr,"mkdir -p %s/%s/toyMC/%s/rev",outPath.c_str(),runName.c_str(),runName.c_str());
      system(tempStr);
      sprintf(tempStr,"%s/%s/toyMC/%s/rev/%s_%i.dat",outPath.c_str(),runName.c_str(),runName.c_str(),runName.c_str(),fEvent);
    }
    std::ofstream toyMCf(tempStr);
    toyMCf << "# x0 y0 zi zf ux uy uz E deltaE deltaX\n";

    double Epart = QMAX;
    for(int i = 0; i < (int) effHits.size()-1; i++) {
      double x = effHits[i]->XYZ().X(); double x1 = effHits[i+1]->XYZ().X();
      double y = effHits[i]->XYZ().Y(); double y1 = effHits[i+1]->XYZ().Y();
      double zi = effHits[i]->XYZ().Z();
      double zf = effHits[i+1]->XYZ().Z();

      double ux = (x1 - x);
      double uy = (y1 - y);
      double uz = (zf - zi);
      double deltaX = sqrt(ux*ux + uy*uy + uz*uz);
      ux /= deltaX; uy /= deltaX; uz /= deltaX;

      double deltaE = effHits[i]->Edep();

      toyMCf << x << " " << y << " " << zi << " " << zf << " " << ux << " " << uy << " " << uz << " " 
             << Epart << " " << deltaE << " " << deltaX << "\n";

      Epart -= deltaE; 
    }
    toyMCf.close();
    
    // Open the output file for this event.
    sprintf(tempStr,"mkdir -p %s/%s/trk",outPath.c_str(),runName.c_str());
    system(tempStr);

    if(forwardFit)
      sprintf(tempStr,"%s/%s/trk/set_f%i.dat",outPath.c_str(),runName.c_str(),fEvent);
    else
      sprintf(tempStr,"%s/%s/trk/set_r%i.dat",outPath.c_str(),runName.c_str(),fEvent);
    std::ofstream outf(tempStr);
    outf << "# node xM yM zM xP yP zP chi2P xF yF zF chi2F\n";

    const std::vector<Node*> tnodes =  trj->nodes();

    auto inode =0;
    for (auto node: tnodes)
    {

      // Declare variables to record position and chi2 information in a file.
      double xM = -1.e9, yM = -1.e9, zM = -1.e9;
      double xP = -1.e9, yP = -1.e9, zP = -1.e9, chi2P = -1.e9;
      double xF = -1.e9, yF = -1.e9, zF = -1.e9, chi2F = -1.e9;

      // Get the measured values from the effective hits.
      if(inode < (int) effHits.size()) {
        xM = effHits[inode]->XYZ().X();
        yM = effHits[inode]->XYZ().Y();
        zM = effHits[inode]->XYZ().Z();
      }
      else {
        std::cout << "WARNING: attempting to access measurement beyond size of eff. hits.";
      }

      // get the state X and C

      State state = node->state();
      //klog << log4cpp::Priority::DEBUG << " *** AVAILABLE KEYS INCLUDE: " << state.hvmap() ;
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
 
      xP = x_P[0]; yP = x_P[1]; zP = x_P[2]; chi2P = tchi2;

      klog << log4cpp::Priority::DEBUG<< " chi2 = " << tchi2 ;

      // Fill the histograms for position and chi2.
      fH1_predictedChi2->Fill(inode,tchi2);
      //fH2_predictedPosition->Fill(x_P[0],x_P[1]);

      // Retrieve filtered residual if present for this node.
      if(state.hvmap().has_key(RP::filtered)) {

        const HyperVector hv_F = state.hv(RP::filtered);  //filtered
        const EVector x_F = hv_F.vector();  // vector state
        const EMatrix C_F = hv_F.matrix();  //cov matrix

        HyperVector resHV_F = node->residuals().hv(RP::filtered); 
        const EVector r_F = resHV_F.vector();
        const EMatrix R_F = resHV_F.matrix();
        double tchi2_F = resHV_F.chi2();

        klog << log4cpp::Priority::DEBUG << " For node " << inode;
        klog << log4cpp::Priority::DEBUG 
        << " X_F[0] " << x_F[0]
        << " X_F[1] " << x_F[1]
        << " X_F[2] " << x_F[2]
        << " X_F[3] " << x_F[3]
        << " X_F[4] " << x_F[4];

        xF = x_F[0]; yF = x_F[1]; zF = x_F[2]; chi2F = tchi2_F;

        klog << log4cpp::Priority::DEBUG<< " chi2 = " << tchi2_F ;

        // Compute the pulls.
        //std::cout << "Pulling with hv_F = " << hv_F << " and pos = " << node->measurement().position_hv() << std::endl;
        //HyperVector residual = hv_F - node->measurement().position_hv();
        EVector pull = resHV_F.pull();
      
        // Fill the histograms. 
        if(tchi2_F < 20)
          fH1_xpull->Fill(pull[0]);   

        fH1_filteredChi2->Fill(inode,tchi2_F);
        //fH2_filteredPosition->Fill(x_F[0],x_F[1]);     
      }

      // Output the information to file.
      outf << inode << " "
           << xM << " " << yM << " " << zM << " " 
           << xP << " " << yP << " " << zP << " " << chi2P << " "
           << xF << " " << yF << " " << zF << " " << chi2F << "\n";

      inode++;
    }

    // Close the output file.
    outf.close();

    fEvent++;

    //getline(cin, input);
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
