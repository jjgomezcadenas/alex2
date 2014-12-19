// A simple class to define a Kalman Filter Setup
// JJ, April 2014

#include <alex/KFSvc.h>
#include <alex/ISvc.h>
#include <iostream>
#include <fstream>
#include <alex/KFSetup.h>
#include <alex/LogUtil.h>
#include <alex/Hit.h>

#include <CLHEP/Units/SystemOfUnits.h> 

using std::string;
using std::endl; 
using std::cout;
using std::cin;
using std::vector;
using std::pair; 

using namespace Recpack;
using namespace CLHEP;

/*
A Kalman Fitter manager is selected by the KFSvc. The fit model can be
either sline (straight line) or helix. The geometrical representation is hard-wired
assuming planes perpendicular to z axis. This, in turn, hardwires the fittin representation
to RP::slopes_z (for sline) or RP::slopes_curv_z (for helix) 
*/

namespace alex {


//--------------------------------------------------------------------
  void KFSvcManager::Init() 
//--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::DEBUG << 
    "In KFSvcManagerSvc::Init() " ;
    
    fModel = Model();
    fDim = ModelDim();
    fKFRep = KFRep();
    
    klog << log4cpp::Priority::INFO 
        << "--: model = " << fModel 
        << " dimension = " << fDim
        << " representation = " << fKFRep;
    
    klog << log4cpp::Priority::DEBUG
        << "Select the model, fitter and rep";

    // select the model
    fRPMan.model_svc().select_model(fModel);

    // select the Kalman Filter as fitting algorithm
    fRPMan.fitting_svc().select_fitter(RP::kalman);

    //select the default representation
    //RP::rep().set_default_rep_name(RP::pos_dir_curv);

     // default fitting representation (should be selected before trying to retrieve a fitter)
    fRPMan.fitting_svc().set_fitting_representation(fKFRep);

    klog << log4cpp::Priority::DEBUG
        << "Select the max chi2, max number outliers & extrap failures";

    // set the maximum local chi2 allowed
    fRPMan.fitting_svc().retrieve_fitter<KalmanFitter>(RP::kalman,fModel).
    set_max_local_chi2ndf(KFSetup::Instance().MaxChi2()); 

    // set the maximum number of outliers allowed 
    fRPMan.fitting_svc().retrieve_fitter<KalmanFitter>(RP::kalman,fModel).
    set_number_allowed_outliers(KFSetup::Instance().MaxOutliers());

    // set the maximum number of consecutive extrapolation failures when predicting in the Kalman Filter
    fRPMan.fitting_svc().retrieve_fitter<KalmanFitter>(RP::kalman,fModel).
    set_max_consecutive_extrap_failures(KFSetup::Instance().MaxExtrapFailures());

    klog << log4cpp::Priority::DEBUG
        << "Init geometrical limits";

    // initialize geometrical limits
    fRPMan.geometry_svc().set_zero_length(1e-3 * mm);
    fRPMan.geometry_svc().set_infinite_length(1e12 * mm);

    // minimum distance between nodes to check the node ordering. 
    // Specially  for P0D which can have clusters at the same layer with different z positions
    fRPMan.matching_svc().set_min_distance_node_ordering(KFSetup::Instance().MinDistanceNodeOrdering());

    // minimum number of nodes to be check for good ordering
    fRPMan.matching_svc().set_min_good_node_ordering(KFSetup::Instance().MinGoodNodeOrdering());

    klog << log4cpp::Priority::DEBUG
        << "Enable MS, disable energy loss fluctuations, enable energy loss correction";

    // enable multiple scattering by default
    fRPMan.model_svc().enable_noiser(fModel, RP::ms, true);

    // enable energy loss fluctuations by default
    fRPMan.model_svc().enable_noiser(fModel, RP::eloss, false);

    // enable electron energy loss fluctuations (bremsstrahlung) by default
    fRPMan.model_svc().enable_noiser(fModel, RP::electron_eloss, false);

    // enable electron energy loss correction (bremsstrahlung) by default
    fRPMan.model_svc().enable_correction(fModel, RP::brem_eloss, false);

    // enable energy loss correction by default
    fRPMan.model_svc().enable_correction(fModel, RP::eloss, false);

    // By default no preselected length sign is used when intersecting a surface
    fRPMan.model_svc().model(fModel).intersector().set_length_sign(0);

    // The geometry is not initialized for this manager
    fRPMan.geometry_svc().set_status("ready",false);

    klog << log4cpp::Priority::DEBUG
        << "Initialize Manager Geometry";

    InitializeManagerGeometry();

    klog << log4cpp::Priority::DEBUG
         << "Set up the surface maker";

    // Set up the surface maker
    fRPMan.matching_svc().add_surface_maker("global",new KFSurfaceMaker(&(fRPMan.matching_svc()), &(fRPMan.geometry_svc().setup())));
    fRPMan.matching_svc().select_surface_maker("global");

    klog << log4cpp::Priority::DEBUG
        << "Set Verbosity();";

    SetVerbosity();

  }
//--------------------------------------------------------------------
  void KFSvcManager::InitializeManagerGeometry()
//--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    
    klog << log4cpp::Priority::DEBUG
        << "Initialize Manager Geometry";

     // box dimensions
    double S = 100*m;

    // the axes for the definition of volumes and surfaces
    
    xaxis = EVector(3,0);
    yaxis = EVector(3,0);
    zaxis = EVector(3,0);
    
    xaxis[0] =1.; yaxis[1]=1.; zaxis[2] = 1.; 

    fBField = EVector(3,0);
    fBField[2]=KFSetup::Instance().B();

    klog << log4cpp::Priority::DEBUG
        << "BField -->" << fBField[2]/tesla;
       
    // 1. Give a name to the setup
    fSetup.set_name("main");

    // 2.  Create a mandatory mother volume
    klog << log4cpp::Priority::DEBUG
        << "Create mother volume";

    EVector pos(3,0);
    pos[0]=0.; pos[1]=0; pos[2]=0;  // at (x,y,z)=(0,0,0)
    Volume* mother = new Box(pos,xaxis,yaxis,S,S,S);
    fSetup.set_mother(mother);

    // 3. Define its properties

    //klog << log4cpp::Priority::DEBUG
    //    << "Use surfaces normal to Z";
    //fSetup.set_volume_property("mother",RP::SurfNormal,zaxis); //DATA MEMBER!!!
    //fSetup.set_volume_property("mother",RP::SurfNormal,xaxis);

    if (KFSetup::Instance().Model() =="helix") 
      fSetup.set_volume_property("mother",RP::BField,fBField);

    if (KFSetup::Instance().X0() != 0) 
    {
      fX0 = KFSetup::Instance().X0(); // Data member, Recpack weird stuff   
      fSetup.set_volume_property("mother",RP::X0,fX0); 
    }

    if (KFSetup::Instance().DEDX() != 0.) 
    {
      fDedx = KFSetup::Instance().DEDX(); 
      fSetup.set_volume_property("mother",RP::de_dx,fDedx);  
    }


    // 4. add the setup to the geometry service
    klog << log4cpp::Priority::DEBUG
        << "Add the setup to the geometry service";

    fRPMan.geometry_svc().add_setup("main",fSetup);

    // 5. select the setup to be used by the geometry service
    fRPMan.geometry_svc().select_setup("main");

     // navigation strategy (the navigator needs the setup)
    klog << log4cpp::Priority::DEBUG
        << "Navigation strategy";

    fRPMan.navigation_svc().navigator(fModel).set_unique_surface(true);

    // The geometry is initialized for this manager
    fRPMan.geometry_svc().set_status("ready",true);

    klog << log4cpp::Priority::INFO 
        << "Manager Geometry: " << fSetup;
  
  }
//--------------------------------------------------------------------
  void KFSvcManager::SetVerbosity() 
//--------------------------------------------------------------------
  {
    //{"MUTE","ERROR","WARNING","INFO","NORMAL","DETAILED","VERBOSE","VVERBOSE"};

    log4cpp::Category& klog = log4cpp::Category::getRoot();

    vector<string> vb =KFSetup::Instance().FiNaMoMaVerbosity();  //FitterNavigationModelMatching
    Messenger::Level level[4];

    for(int i=0; i <4; ++i)
    {
      level[i] =Messenger::str(vb.at(i));

      klog << log4cpp::Priority::DEBUG 
        << " In SetRPVerbosity: \n  " 
        << " level: " << i << " verbosity = " << vb.at(i) <<  "\n";
    }

    klog << log4cpp::Priority::DEBUG 
        << " In SetRPVerbosity: \n  " 
        << " level0: " << level[0] << "\n"
        << " level1: " << level[1] << "\n"
        << " level2: " << level[2] << "\n"
        << " level3: " << level[3] << "\n";
        

    Messenger::Level level0 = level[0];
    Messenger::Level level1 = level[1];
    Messenger::Level level2 = level[2];
    Messenger::Level level3 = level[3];
    Messenger::Level level4 = level[3];

    const std::string& modelname = fRPMan.model_svc().model_name();

    klog << log4cpp::Priority::DEBUG 
        << " model name = " << modelname ;
        
    klog << log4cpp::Priority::DEBUG 
        << " set verbosity for fitting: level0 = " << level0 ;

    fRPMan.fitting_svc().fitter(modelname).set_verbosity(level0);

    klog << log4cpp::Priority::DEBUG 
        << " set verbosity for navigation: level1 = " << level1 ;

    fRPMan.navigation_svc().set_verbosity(level1);
    fRPMan.navigation_svc().navigator(modelname).set_verbosity(level1);
    fRPMan.navigation_svc().inspector(RP::X0).set_verbosity(level1);
    fRPMan.navigation_svc().navigator(modelname).master_inspector().set_verbosity(level1);
    fRPMan.navigation_svc().inspector(RP::BField).set_verbosity(level1);  
    fRPMan.navigation_svc().inspector(RP::eloss).set_verbosity(level1);     
    fRPMan.navigation_svc().inspector(RP::Nel).set_verbosity(level1);     
    fRPMan.navigation_svc().inspector(RP::elossMap).set_verbosity(level1);     
    fRPMan.navigation_svc().inspector(RP::BFieldMap).set_verbosity(level1);     
    fRPMan.model_svc().model(modelname).intersector("plane").set_verbosity(level4);
    fRPMan.model_svc().model(modelname).intersector("numerical").set_verbosity(level1);

    klog << log4cpp::Priority::DEBUG 
        << " set verbosity for model: level2 = " << level2;

    fRPMan.model_svc().model(modelname).equation().set_verbosity(level2);
    fRPMan.model_svc().model(modelname).propagator().set_verbosity(level2);
    fRPMan.model_svc().model(modelname).noiser().set_verbosity(level2);
    fRPMan.model_svc().model(modelname).tool("correction/eloss").set_verbosity(level2);
    fRPMan.model_svc().model(modelname).tool("correction/brem_eloss").set_verbosity(level2);
    fRPMan.model_svc().model(modelname).tool("noiser/brem_eloss").set_verbosity(level2);
    fRPMan.model_svc().model(modelname).tool("noiser/ms").set_verbosity(level2);
    fRPMan.model_svc().model(modelname).tool("noiser/eloss").set_verbosity(level2);
    fRPMan.model_svc().model(modelname).tool("intersector/numerical").set_verbosity(level2);

    klog << log4cpp::Priority::DEBUG 
        << " set verbosity for matching: level3 = " << level3;

    fRPMan.conversion_svc().projector().set_verbosity(level2);
    fRPMan.conversion_svc().representation(RP::slopes_curv).set_verbosity(level2);

    fRPMan.matching_svc().set_verbosity(level3);
  //fRPMan.matching_svc().surface_maker("global").set_verbosity(level3);

    RayTool::msg().set_verbosity(level4);

  }

//--------------------------------------------------------------------
  RP::Trajectory* KFSvcManager::CreateTrajectory(std::vector<const Hit* > hits, 
                                                std::vector<double> hitErrors) 
//--------------------------------------------------------------------
  {
    
    log4cpp::Category& klog = log4cpp::Category::getRoot();

    klog << log4cpp::Priority::DEBUG << " In CreateTrajectory --> " ;
  
    klog << log4cpp::Priority::DEBUG << "Cov Matrix --> " ;
    EVector m(3,0);
    fCov = EMatrix(3,3,0);
    
    for (size_t i=0; i<3; ++i)
    {
      fCov[i][i] = hitErrors[i];
    }

    klog << log4cpp::Priority::DEBUG << "Cov Matrix --> " << fCov;

    // destroy measurements in fMeas container
    //stc_tools::destroy(fMeas);

    RP::measurement_vector fMeas;

    auto size = hits.size();

    klog << log4cpp::Priority::DEBUG << " size of hits --> " << size;
    klog << log4cpp::Priority::DEBUG << " fill measurement vector --> " << size;

    for(auto hit : hits)
    {
      m[0]=hit->XYZ().X();
      m[1]=hit->XYZ().Y();
      m[2]=hit->XYZ().Z();
 
      klog << log4cpp::Priority::DEBUG << "+++hit x = " << m[0] 
      << " y = " << m[1] << " z = " << m[2];
      klog << log4cpp::Priority::DEBUG << "Create a measurement " ;

      Measurement* meas = new Measurement(); 
      string type=RP::xyz;
   
      meas->set_name(type);                         // the measurement type
      meas->set_hv(HyperVector(m,fCov,type));  // the HyperVector 
      meas->set_deposited_energy(hit->Edep()); // the deposited energy
      meas->set_name(RP::setup_volume,"mother");          // the volume
        
      // the  position of the measurement
      klog << log4cpp::Priority::DEBUG << "Create measurement plane " ;

      meas->set_position_hv(HyperVector(m,EMatrix(3,3,0),RP::xyz));

      //    Add the measurement to the vector
      fMeas.push_back(meas);
    }
    
    // add the measurements to the trajectory
    klog << log4cpp::Priority::INFO << " Measurement vector created "  ;
    
    Trajectory* trj = new Trajectory();
    trj->add_constituents(fMeas);  

    klog << log4cpp::Priority::INFO << " Trajectory --> " << *trj ;

    return trj;
  }

//--------------------------------------------------------------------
  RP::State* KFSvcManager::SeedState(std::vector<double> v0, std::vector<double> p0) 
//--------------------------------------------------------------------
  {
    
    State* state = new State();
    // v0 is a guess of the position
    //p0 is a guess of the momentum

    log4cpp::Category& klog = log4cpp::Category::getRoot();
    
    EVector v(fDim,0);
    EMatrix C(fDim,fDim,0);   

    v[0] = v0[0];   // a guess of the x position
    v[1] = v0[1];   // a guess of the y position
    v[2] = v0[2];  // a guess of the z position  
    v[3] = p0[0]/p0[2];  // a guess of dx/dz 
    v[4] = p0[1]/p0[2];  // a guess of dy/dz 

    klog << log4cpp::Priority::INFO << " Initial state --> " 
          << v;

    double pmag = sqrt(p0[0]*p0[0] + p0[1]*p0[1] + p0[2]*p0[2]);
    double qoverp = -1./pmag;

    if (KFSetup::Instance().Model() =="helix" && KFSetup::Instance().FitMomentum())
    {
      v[5]=qoverp;  // assume forward fit!
      C[5][5]= 0.1;  // not a large error like the others, momentum "known"
      
      // Set a qoverp HyperVector for use only in multiple scattering calculations.
      klog << log4cpp::Priority::INFO << " Helix model, set qoverp --> " << qoverp ;
      state->set_hv("MSqoverp",HyperVector(qoverp,0));
      state->set_hv("MSdeltaqoverp",HyperVector(0.,0));
    }
    else
    {
      // For the Straight line model q/p is a fix parameter use for MS computation 
    
      klog << log4cpp::Priority::INFO << " Sline model, set qoverp --> " << qoverp ;
      state->set_hv(RP::qoverp,HyperVector(qoverp,0));
      std::cout << "qoverp hypervector is " << state->hv(RP::qoverp).vector()[0] << std::endl;
    }

    // diagonal covariance matrix
    C[0][0] = C[1][1] = 0.1*cm;
    C[2][2] = EGeo::zero_cov()/2.; // no error on Z since planes are perpendicular to z
    //C[0][0] = C[1][1] = C[2][2] = 100.*cm;
    C[3][3] = C[4][4] = 1.;

    klog << log4cpp::Priority::INFO << " Cov matrix --> " << C ;

    // Set the main HyperVector
    klog << log4cpp::Priority::INFO << " Set the main HyperVector --> "  ;
    state->set_hv(HyperVector(v,C));

    // Set the sense HyperVector (1=increasing z)
    double sense=p0[2]/fabs(p0[2]); 

    klog << log4cpp::Priority::INFO << " Set the sense HyperVector --> " << sense ;
    state->set_hv(RP::sense,HyperVector(sense,0));

    // Set the model name
    klog << log4cpp::Priority::INFO << " Set the model name --> " << fModel ;
    state->set_name(fModel);

    // Set the representation (x,y,z dx/dz, dy/dz, q/p)
    klog << log4cpp::Priority::INFO << " Set the representation --> "  ;
    if (KFSetup::Instance().Model() =="helix" && KFSetup::Instance().FitMomentum())
      state->set_name(RP::representation,RP::slopes_curv_z);
    else
      state->set_name(RP::representation,RP::slopes_z);

    // Set the particle type to electron.
    state->set_name(RP::PID, "Electron");

    return state;
  }

// Returns the number of nodes on which the fit broke.
//--------------------------------------------------------------------
  int KFSvcManager::FitTrajectory(RP::Trajectory& traj, RP::State& seed, std::vector<int>& fit_breaks) 
//--------------------------------------------------------------------
  { 
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::INFO << "In KFSvcManager::Fitting -- " ;

    // Clear the list of fit breaks.
    fit_breaks.clear();

    // ---------------------------
    // Set up for the initial fit.
    // ---------------------------
    bool status = false;
    int last_fitted = 0;
    int tot_nodes = traj.nodes().size();

    // Perform the initial fit.
    status = fRPMan.fitting_svc().fit(seed,traj);
    last_fitted = traj.last_fitted_node();
    //std::cout << "-- [FIT 0] returned status " << status << std::endl;

    // If the entire fit succeeded, return.
    if(last_fitted == tot_nodes-1) return 0; // status;

    // Add the last fit node to the list of fit breaks.
    fit_breaks.push_back(last_fitted);

    // If the fit did not succeed completely, fit again using a new trajectory constructed from
    //  the remaining nodes that were not fit.
   
    // Keep track of the energy lost.
    double elost = 0.;
    
    // Copy the trajectory and clear the original.
    RP::Trajectory fitTraj(traj);
    traj.reset();

    // Loop until entire track is fit or maximum number of breaks has been reached.
    int num_breaks = 1;
    bool all_fit = false;
    while(!all_fit && num_breaks < 40) {

        // Add all fit nodes to the final trajectory and keep track of the energy lost.
        const std::vector<Node*> tnodes =  fitTraj.nodes();
        for(int n = 0; n < fitTraj.last_fitted_node(); n++) {
            traj.add_node(*tnodes[n]);
            elost += tnodes[n]->measurement().deposited_energy();
        }
        //std::cout << "-- [FIT " << num_breaks-1 << "]: " << fitTraj.last_fitted_node()+1 << " nodes" << std::endl;

        // Create a new trajectory from the nodes that were not fit.
        std::vector<Node*> tofit_nodes;
        for(int n = fitTraj.last_fitted_node(); n < (int) tnodes.size(); n++) {
            Node * nnode = new Node(*tnodes[n]);
            tofit_nodes.push_back(nnode);
        }
        fitTraj.reset();
        fitTraj.add_nodes(tofit_nodes);

        // Only attempt to fit again if more than 2 nodes are left.
        if(tofit_nodes.size() < 2) {
            status = true;
            all_fit = true;
        }
        else {
       
            //std::cout << "Creating new seed state..." << std::endl;
 
            // Create a new seed state.
            std::vector<double> v0;
            v0.push_back(tofit_nodes[0]->measurement().position()[0]);
            v0.push_back(tofit_nodes[0]->measurement().position()[1]);
            v0.push_back(tofit_nodes[0]->measurement().position()[2]);

            std::vector<double> p0;
            std::vector<double> xlist; std::vector<double> ylist; std::vector<double> zlist;
            xlist.push_back(tofit_nodes[0]->measurement().position()[0]); 
            ylist.push_back(tofit_nodes[0]->measurement().position()[1]); 
            zlist.push_back(tofit_nodes[0]->measurement().position()[2]);

            xlist.push_back(tofit_nodes[1]->measurement().position()[0]);
            ylist.push_back(tofit_nodes[1]->measurement().position()[1]); 
            zlist.push_back(tofit_nodes[1]->measurement().position()[2]);        
            ISvc::Instance().GuessInitialMomentum(QMAX+MELEC-elost, p0, xlist, ylist, zlist);
 
            State * fitSeed = SeedState(v0,p0);

            // Perform an additional fit.
            status = fRPMan.fitting_svc().fit(*fitSeed,fitTraj);

            if(fitTraj.last_fitted_node() == ((int)fitTraj.nodes().size())-1) all_fit = true;
            else fit_breaks.push_back(fitTraj.last_fitted_node() + traj.nodes().size());
        }

        num_breaks++;
    }

    // Add the remaining fit nodes to the main track.
    const std::vector<Node*> fnodes =  fitTraj.nodes();
    //for(int n = 0; n < fitTraj.last_fitted_node(); n++) {
    for(int n = 0; n < (int) fitTraj.nodes().size(); n++) {
        traj.add_node(*fnodes[n]);
    }

    // Ensure we have fit all nodes.
    if(!all_fit) status = false;

    return num_breaks;
      
  }
//--------------------------------------------------------------------
  int KFSvcManager::ModelDim() const 
//--------------------------------------------------------------------
  {
    int dim;
    if (KFSetup::Instance().Model() =="helix" && KFSetup::Instance().FitMomentum())
    {
      dim = 6;                     // dimension of the state vector for this model
    } 
    else // if (KFSetup::Instance().Model() =="sline")
    {
      dim = 5; 
    }
    //else
    //{
    //  std::cout << "ERROR: model unknown" << std::endl;
    //  exit (EXIT_FAILURE);
    //}
    return dim;
  }
//--------------------------------------------------------------------
  string KFSvcManager::Model() const 
//--------------------------------------------------------------------
  {
    string model;
    if (KFSetup::Instance().Model() =="helix")
    {
      model = RP::particle_helix;         // model used for reconstruction                   
    } 
    else if (KFSetup::Instance().Model() =="sline")
    {
      model = RP::particle_sline;         // model used for reconstruction
    }
    else
    {
      std::cout << "ERROR: model unknown" << std::endl;
      exit (EXIT_FAILURE);
    }
    return model;
  }
  //--------------------------------------------------------------------
  string KFSvcManager::KFRep() const 
 //--------------------------------------------------------------------
  {
    string kfrep;

    if (KFSetup::Instance().Model() =="helix" && KFSetup::Instance().FitMomentum())
    {
      kfrep =  RP::slopes_curv;        // (x,y,z,ux,uy,q/p)  // RP::slopes_curv_z
    } 
    else //if (KFSetup::Instance().Model() =="sline")
    {
      kfrep = RP::slopes;
      //kfrep =  RP::slopes_z ;        // (x,y,z,ux,uy)
    }
    //else
    //{
    //  std::cout << "ERROR: model unknown" << std::endl;
    //  exit (EXIT_FAILURE);
    //}
    
    return kfrep;
  }

  //-------------------------------------------------------------------------------------------------------------------------------
  void KFSvcManager::OutputFit(const char * fname, const std::vector<Node*> & nodes, const std::vector<int> & fit_breaks) const {
  //-------------------------------------------------------------------------------------------------------------------------------

    std::cout << "Outputting fit to file " << fname << std::endl;

    std::ofstream outf(fname);
    outf << "# node xM yM zM xP yP zP chi2P xF yF zF chi2F qoverp qoverpfit deltaE brk sense\n";
          
    auto inode = 0;
    for (auto node: nodes)  {

      // Declare variables to record position and chi2 information in a file.
      double xM = -1.e9, yM = -1.e9, zM = -1.e9;
      double xP = -1.e9, yP = -1.e9, zP = -1.e9, chi2P = -1.e9;
      double xF = -1.e9, yF = -1.e9, zF = -1.e9, chi2F = -1.e9;
      double eM = -1.e9;
      double varqoverp = -1.e9;
      double stqoverp = -1.e9;
      int fitsense = -1.e9;
      int brk = 0;

      // Get the measured values.
      xM = node->measurement().position()[0];   // effHits[lseg_i + lseg_inode]->XYZ().X();
      yM = node->measurement().position()[1];   // effHits[lseg_i + lseg_inode]->XYZ().Y();
      zM = node->measurement().position()[2];   // effHits[lseg_i + lseg_inode]->XYZ().Z();
      eM = node->measurement().deposited_energy(); // effHits[lseg_i + lseg_inode]->Edep();

      // Determine if the fit was broken just after fitting this node.
      for(unsigned int ibrk = 0; ibrk < fit_breaks.size(); ibrk++) {
        if(fit_breaks[ibrk] == inode) brk = 1;
      }
        
      if(node->status(RP::fitted)) {

        State state = node->state();
        const HyperVector hv_P = state.hv(RP::predicted);  //predicted 
        const EVector x_P = hv_P.vector();  // vector state
        const EMatrix C_P = hv_P.matrix();  //cov matrix

        // get the q/p
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
        HyperVector resHV_P = node->residuals().hv(RP::predicted);
        const EVector r_P = resHV_P.vector();
        const EMatrix R_P = resHV_P.matrix();
        double tchi2 = resHV_P.chi2();
        xP = x_P[0]; yP = x_P[1]; zP = x_P[2]; chi2P = tchi2;

        // Retrieve filtered residual if present for this node.
        if(state.hvmap().has_key(RP::filtered)) {
          
          const HyperVector hv_F = state.hv(RP::filtered);  //filtered
          const EVector x_F = hv_F.vector();  // vector state
          const EMatrix C_F = hv_F.matrix();  //cov matrix
       
          HyperVector resHV_F = node->residuals().hv(RP::filtered);
          const EVector r_F = resHV_F.vector();
          const EMatrix R_F = resHV_F.matrix();
          double tchi2_F = resHV_F.chi2();

          xF = x_F[0]; yF = x_F[1]; zF = x_F[2]; chi2F = tchi2_F;
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
    }
    std::cout << "End of output function." << std::endl;
  }
}
