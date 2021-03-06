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

    KFSetup::Instance().Init(gas, model, Pr, B, fFitMomentum);
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

    // Get the event number for this event.
    fEvent = ISvc::Instance().GetEvtNum();

    // Compute the total energy in the event.
    double etot = 0.;
    IHits truehits = ISvc::Instance().GetTrueHits();
    for(int h = 0; h < (int) truehits.size(); h++) {
      etot += truehits[h].second;
    } 

    //RBeta* rBeta = ISvc::Instance().GetRBeta();
    //std::vector<const Hit*> effHits = rBeta->GetEffHits();

    bool forwardFit = true;
    if(fForwardFit == 0) forwardFit = false;

    // Fit assuming the initial kinetic energy is QMAX.
    ISvc::Instance().CreateKFObjects(QMAX+MELEC,errXY,true,forwardFit);

    std::vector<const Hit*> effHits = ISvc::Instance().GetKFHits();
    std::vector<double> hitErrors = ISvc::Instance().GetKFMErrors();
    std::vector<double> v0 = ISvc::Instance().GetKFv0();
    std::vector<double> p0 = ISvc::Instance().GetKFp0();

    // Add the obtained effective hits to a histogram.
    double eeff = 0.;
    for(int i = 0; i < (int) effHits.size(); i++) {
      eeff += effHits[i]->Edep();
      //double x = effHits[i]->XYZ().X();
      //double y = effHits[i]->XYZ().Y();
      //fH2_xyMeasured->Fill(x,y,i);
    }

    fH1_energy->Fill(eeff);
    fH1_ediff->Fill(etot-eeff);
    fH1_momentum->Fill(sqrt((etot + 0.511)*(etot + 0.511) - 0.511*0.511));

    string input ;

    klog << log4cpp::Priority::DEBUG << " Creating trajectory" ;
    //getline(cin, input);
   
    Trajectory* trj = KFSvc::Instance().CreateTrajectory(effHits, hitErrors);
                                      
    klog << log4cpp::Priority::DEBUG << " Creating seed" ;
    //getline(cin, input);

    State* seed =KFSvc::Instance().SeedState(v0, p0) ;

   
    klog << log4cpp::Priority::DEBUG << " Fitting trajectory" ;
    //getline(cin, input);

    //bool status = KFSvc::Instance().FitTrajectory(*trj, *seed);
    std::vector<int> fit_breaks;
    int nbreaks = KFSvc::Instance().FitTrajectory(*trj, *seed, fit_breaks);

    //klog << log4cpp::Priority::DEBUG << " Fitting trajectory, result " << status;
    klog << log4cpp::Priority::DEBUG << " Fitting trajectory, result " << nbreaks << " fits";
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
      sprintf(tempStr,"%s/%s/trk/%s_f%i.dat",outPath.c_str(),runName.c_str(),runName.c_str(),fEvent);
    else
      sprintf(tempStr,"%s/%s/trk/%s_r%i.dat",outPath.c_str(),runName.c_str(),runName.c_str(),fEvent);
    std::ofstream outf(tempStr);
    outf << "# node xM yM zM xP yP zP chi2P xF yF zF chi2F qoverp qoverpfit deltaE brk sense\n";

    const std::vector<Node*> tnodes =  trj->nodes();

    auto inode =0;
    for (auto node: tnodes)
    {

      // Declare variables to record position and chi2 information in a file.
      double xM = -1.e9, yM = -1.e9, zM = -1.e9;
      double xP = -1.e9, yP = -1.e9, zP = -1.e9, chi2P = -1.e9;
      double xF = -1.e9, yF = -1.e9, zF = -1.e9, chi2F = -1.e9;
      double eM = -1.e9;
      double Estate = 0.;
      double varqoverp = -1.e9;
      double stqoverp = -1.e9;
      int fitsense = -1.e9;
      int brk = 0;

      // Get the measured values from the effective hits.
      if(inode < (int) effHits.size()) {
        xM = effHits[inode]->XYZ().X();
        yM = effHits[inode]->XYZ().Y();
        zM = effHits[inode]->XYZ().Z();
        eM = effHits[inode]->Edep();
      }
      else {
        std::cout << "WARNING: attempting to access measurement beyond size of eff. hits.";
      }

      // Determine if the fit was broken just after fitting this node.
      for(unsigned int ibrk = 0; ibrk < fit_breaks.size(); ibrk++) {
       if(fit_breaks[ibrk] == inode) brk = 1;  
      }

      if(node->status(RP::fitted)) {

        State state = node->state();
        //klog << log4cpp::Priority::DEBUG << " *** AVAILABLE KEYS INCLUDE: " << state.hvmap() ;
        const HyperVector hv_P = state.hv(RP::predicted);  //predicted 
        const EVector x_P = hv_P.vector();  // vector state
        const EMatrix C_P = hv_P.matrix();  //cov matrix

        // get the q/p
        int index_qoverp;
        bool ok = Recpack::rep_default().index(RP::qoverp, index_qoverp);
        if(ok) std::cout << "Fitted q/p = " << index_qoverp <<  std::endl;
        if(state.hvmap().has_key("qoverp")) {
          varqoverp = state.hv("qoverp").vector()[0];
        }
        else if(state.hvmap().has_key("MSqoverp")) {
          varqoverp = state.hv("MSqoverp").vector()[0];
          stqoverp = x_P[5];
        }
        else {
          std::cout << "WARNING: no qoverp or MSqoverp found.";
        }

        // get the sense
        fitsense = state.hv(RP::sense).vector()[0];

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
      }

      // Output the information to file.
      outf << inode << " "
           << xM << " " << yM << " " << zM << " " 
           << xP << " " << yP << " " << zP << " " << chi2P << " "
           << xF << " " << yF << " " << zF << " " << chi2F << " "
           << varqoverp << " " << stqoverp << " " << eM << " "
           << brk << " " << fitsense << "\n";

      inode++;
    }

    // Close the output file.
    outf.close();

    // -------------------------------------------------------------------------------------
    // Perform forward and reverse fits for the longest trajectory if this is a forward fit.
    // -------------------------------------------------------------------------------------
    if(forwardFit && fit_breaks.size() > 0) {

      // First find the beginning and end of the largest segment.
      int lseg_i = 0, lseg_f = -1;
      int lseg_max = -1;
      for(unsigned int si = 0; si < fit_breaks.size(); si++) {

        int lseg_size = -1;
     
        // Treat the case of the first break. 
        if(si == 0) {
          lseg_size = fit_breaks[0] + 1;
          lseg_max = lseg_size;
          lseg_i = 0; lseg_f = fit_breaks[0]; 
        }
        else {
  
          lseg_size = fit_breaks[si] - fit_breaks[si-1];

          // Update the indices if we have found a new largest segment.
          if(lseg_size > lseg_max) {
            lseg_i = fit_breaks[si-1]; lseg_f = fit_breaks[si];
            lseg_max = lseg_size;
          }
        }
      }

      std::cout << "Found longest segment with lseg_i = " << lseg_i << " and lseg_f = " << lseg_f << std::endl;

      // Set up variables used to construct and fit longest segment.
      std::vector<Node*> lseg_fnodes; std::vector<Node*> lseg_rnodes;
      std::vector<double> lseg_v0;
      std::vector<double> lseg_p0;
      std::vector<double> lseg_xlist; std::vector<double> lseg_ylist; std::vector<double> lseg_zlist;
 
      // Build the segment in the forward direction.
      double elostf = 0.;
      for(int n = lseg_i; n < lseg_f; n++) {
        if(n >= lseg_i) {
          Node * nnode = new Node(*tnodes[n]);
          lseg_fnodes.push_back(nnode);
        }
        else {
          elostf += tnodes[n]->measurement().deposited_energy();
        }
      }

      // Build the segment in the reverse direction.
      double elostr = 0.;
      for(int n = (int) (tnodes.size()-1); n > lseg_i; n--) {
        std::cout << "Looking to add nodes with n = " << n << std::endl;
        if(n < lseg_f) {
          Node * nnode = new Node(*tnodes[n]);
          lseg_rnodes.push_back(nnode);
          std::cout << "Added node " << n << std::endl;
        }
        else {
          elostr += tnodes[n]->measurement().deposited_energy();
        }
      }

      // Forward fit to longest segment --------------------------------------------------------------------------

      trj->reset();
      trj->add_nodes(lseg_fnodes);

      std::cout << "Built forward segment with " << lseg_fnodes.size() << " nodes; elost = " << elostf << std::endl;

      // Create a seed state for the segment.
      lseg_v0.push_back(lseg_fnodes[0]->measurement().position()[0]);
      lseg_v0.push_back(lseg_fnodes[0]->measurement().position()[1]);
      lseg_v0.push_back(lseg_fnodes[0]->measurement().position()[2]);

      lseg_xlist.push_back(lseg_fnodes[0]->measurement().position()[0]);
      lseg_ylist.push_back(lseg_fnodes[0]->measurement().position()[1]);
      lseg_zlist.push_back(lseg_fnodes[0]->measurement().position()[2]);

      lseg_xlist.push_back(lseg_fnodes[1]->measurement().position()[0]);
      lseg_ylist.push_back(lseg_fnodes[1]->measurement().position()[1]);
      lseg_zlist.push_back(lseg_fnodes[1]->measurement().position()[2]);
      ISvc::Instance().GuessInitialMomentum(QMAX+MELEC-elostf, lseg_p0, lseg_xlist, lseg_ylist, lseg_zlist);
 
      State * lseg_seedf = KFSvc::Instance().SeedState(lseg_v0,lseg_p0);

      // Fit the segment.
      fit_breaks.clear();
      nbreaks = KFSvc::Instance().FitTrajectory(*trj, *lseg_seedf, fit_breaks);

      std::cout << "Forward fit completed with " << nbreaks << " breaks" << std::endl;

      // Output the results.
      sprintf(tempStr,"%s/%s/trk/%s_flsegf%i.dat",outPath.c_str(),runName.c_str(),runName.c_str(),fEvent);
      KFSvc::Instance().OutputFit(tempStr, trj->nodes(), fit_breaks);

      // Clear the temporary variables/arrays ----------------------------------------------------------------------
      elostr = 0.;
      for(int n = 0; n < (int) lseg_fnodes.size(); n++) delete lseg_fnodes[n];
      lseg_v0.clear(); lseg_p0.clear();
      lseg_xlist.clear(); lseg_ylist.clear(); lseg_zlist.clear();

      // Reverse fit to longest segment ----------------------------------------------------------------------------

      std::cout << "Building reverse segment " << std::endl;

      trj->reset();
      trj->add_nodes(lseg_rnodes);

      std::cout << "Built reverse segment with " << lseg_rnodes.size() << " nodes; elost = " << elostr << std::endl;

      // Create a seed state for the segment.
      lseg_v0.push_back(lseg_rnodes[0]->measurement().position()[0]);
      lseg_v0.push_back(lseg_rnodes[0]->measurement().position()[1]);
      lseg_v0.push_back(lseg_rnodes[0]->measurement().position()[2]);

      lseg_xlist.push_back(lseg_rnodes[0]->measurement().position()[0]);
      lseg_ylist.push_back(lseg_rnodes[0]->measurement().position()[1]);
      lseg_zlist.push_back(lseg_rnodes[0]->measurement().position()[2]);

      lseg_xlist.push_back(lseg_rnodes[1]->measurement().position()[0]);
      lseg_ylist.push_back(lseg_rnodes[1]->measurement().position()[1]);
      lseg_zlist.push_back(lseg_rnodes[1]->measurement().position()[2]);
      ISvc::Instance().GuessInitialMomentum(QMAX+MELEC-elostr, lseg_p0, lseg_xlist, lseg_ylist, lseg_zlist);
 
      State * lseg_seedr = KFSvc::Instance().SeedState(lseg_v0,lseg_p0);

      // Fit the segment.
      fit_breaks.clear();
      nbreaks = KFSvc::Instance().FitTrajectory(*trj, *lseg_seedr, fit_breaks);
      std::cout << "Reverse fit completed with " << nbreaks << " breaks" << std::endl;

      // Output the results.    
      sprintf(tempStr,"%s/%s/trk/%s_flsegr%i.dat",outPath.c_str(),runName.c_str(),runName.c_str(),fEvent);
      KFSvc::Instance().OutputFit(tempStr, trj->nodes(), fit_breaks);

      for(int n = 0; n < (int) lseg_rnodes.size(); n++) delete lseg_rnodes[n];
      delete lseg_seedf;
      delete lseg_seedr;
    }

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
