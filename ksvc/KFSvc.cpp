// A simple class to define a Kalman Filter Setup
// JJ, April 2014

#include "KFSvcManager.h"
#include <iostream>
#include <alex/KFSetup.h>
#include <alex/LogUtil.h>
#include <alex/ISvc.h>
#include <CLHEP/Units/SystemOfUnits.h>

using std::string;
using std::endl; 
using std::cout;
using std::cin;
using std::vector;
using std::pair; 

using namespace Recpack;
using namespace CLHEP;

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
    
    
    klog << log4cpp::Priority::INFO 
        << "--: model = " << fModel 
        << " dimension = " << fDim;
    
    klog << log4cpp::Priority::DEBUG
        << "Select the model, fitter and rep";

    fRPMan.model_svc().select_model(fModel);
    fRPMan.fitting_svc().select_fitter(RP::kalman);
    RP::rep().set_default_rep_name(RP::pos_dir_curv);
    fRPMan.fitting_svc().set_fitting_representation(RP::slopes_curv_z);

    klog << log4cpp::Priority::DEBUG
        << "Select the max chi2, max number outliers & extrap failures";

    fRPMan.fitting_svc().retrieve_fitter<KalmanFitter>(RP::kalman,fModel).
    set_max_local_chi2ndf(KFSetup::Instance().MaxChi2());  
    fRPMan.fitting_svc().retrieve_fitter<KalmanFitter>(RP::kalman,fModel).
    set_number_allowed_outliers(KFSetup::Instance().MaxOutliers());
    fRPMan.fitting_svc().retrieve_fitter<KalmanFitter>(RP::kalman,fModel).
    set_max_consecutive_extrap_failures(KFSetup::Instance().MaxExtrapFailures());

    klog << log4cpp::Priority::DEBUG
        << "Init geometrical limits";
    fRPMan.geometry_svc().set_zero_length(1e-3 * mm);
    fRPMan.geometry_svc().set_infinite_length(1e12 * mm);

    // minimum distance between nodes to check the node ordering. 
    // Specially  for P0D which can have clusters at the same layer with different z positions
    fRPMan.matching_svc().set_min_distance_node_ordering(fKfs->MinDistanceNodeOrdering());

    // minimum number of nodes to be check for good ordering
    fRPMan.matching_svc().set_min_good_node_ordering(KFSetup::Instance().MinGoodNodeOrdering());

    klog << log4cpp::Priority::DEBUG
        << "Enable MS, disable energy loss fluctuations, enable energy loss correction";
    fRPMan.model_svc().enable_noiser(fModel, RP::ms, true);

    // disable energy loss fluctuations by default
    fRPMan.model_svc().enable_noiser(fModel, RP::eloss, false);

    // disable electron energy loss fluctuations (bremsstrahlung) by default
    fRPMan.model_svc().enable_noiser(fModel, RP::electron_eloss, false);

    // disable electron energy loss correction (bremsstrahlung) by default
    fRPMan.model_svc().enable_correction(fModel, RP::brem_eloss, false);

    // enable energy loss correction by default
    fRPMan.model_svc().enable_correction(fModel, RP::eloss, true);

    // By default no preselected length sign is used when intersecting a surface
    fRPMan.model_svc().model(fModel).intersector().set_length_sign(0);

    // The geometry is not initialized for this manager
    fRPMan.geometry_svc().set_status("ready",false);

    klog << log4cpp::Priority::DEBUG
        << "Create the setup";
    CreateSetup();

   klog << log4cpp::Priority::DEBUG
        << "Setup Created";
    SetVerbosity();

  }
//--------------------------------------------------------------------
  void KFSvcManager::CreateSetup()
//--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    
    klog << log4cpp::Priority::DEBUG
        << "Create CreateKFSetup";

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

    klog << log4cpp::Priority::DEBUG
        << "Use surfaces normal to Z";
    fSetup.set_volume_property("mother",RP::SurfNormal,zaxis); //DATA MEMBER!!!

    if (KFSetup::Instance().Model() =="helix") 
      fSetup.set_volume_property("mother",RP::BField,fBField);

    if (KFSetup::Instance().X0() != 0) 
    {
      double fX0 = KFSetup::Instance().X0(); // Data member, Recpack weird stuff   
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
        << "Recpack Setup: " << fSetup;
  
    klog << log4cpp::Priority::INFO 
        << "Dump completed: " ;
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
    fRPMan.model_svc().model(modelname).intersector("plane").set_verbosity(level1);
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
  RP::Trajectory KFSvcManager::CreateTrajectory() 
//--------------------------------------------------------------------
  {
    
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::INFO << "+++In KFSvcManager::Create Trajectory " ;

    std::vector<const IBeta*> iBetas = ISvc::Instance().GetIBetas()

    for(auto beta : iBetas)
    {
      if (not beta->GetType()=="2e")
        continue;

      
    }

    IHits sphits = isvc.SelectHits("Sparse");
    IHits smhits = isvc.SelectHits("Smeared");
    std::vector<double> sigma = isvc.MeasurementErrors();
    std::vector<const irene::Particle*> vBeta = isvc.Electrons();

    const irene::Particle* beta = vBeta[0];
    auto betaTest = isvc.SelectBetaP(beta);    // it test= true then single electron

    klog << log4cpp::Priority::INFO << "SelectBetaP --> " << betaTest;


    // destroy measurements in _meas container
    stc_tools::destroy(fMeas);
  
    //reset the "true" trajectory
    fSimTraj.reset();
  
    EVector vv(7,0);
    EMatrix M(3,3,1);
    EVector m(3,0);
    EVector smear_m(3,0);
    EMatrix meas_cov(3,3,0);
    
    for (size_t i=0; i<3; ++i)
    {
      meas_cov[i][i] = sigma[i];
    }

    klog << log4cpp::Priority::DEBUG << "Cov Matrix --> " << meas_cov;

    auto size = sphits.size();

    if (size != smhits.size())
    {
      std::cout << "ERROR: sparse and smeared hits do not have the same size" << std::endl;
      std::cout << "sphits.size() =" << sphits.size() << std::endl;
      std::cout << "smhits.size() =" << smhits.size() << std::endl;
      exit (EXIT_FAILURE);
    }

    klog << log4cpp::Priority::DEBUG << " size of hits --> " << size;

    for (size_t i=0; i<size; ++i)
    {
      auto sphit =sphits.at(i);
      auto smhit =smhits.at(i); 
      auto xyztSp = sphit.first;
      auto xyztSm = smhit.first;

      m[0]=xyztSp.X();
      m[1]=xyztSp.Y();
      m[2]=xyztSp.Z();

      smear_m[0]=xyztSm.X();
      smear_m[1]=xyztSm.Y();
      smear_m[2]=xyztSm.Z();
    
      string type=RP::xyz;

      klog << log4cpp::Priority::DEBUG << "+++true hit -> " << m;
      klog << log4cpp::Priority::DEBUG << "+++smear hit -> " << smear_m;
      
      // create a measurement
      klog << log4cpp::Priority::DEBUG << "Create a measurement " ;
      Measurement* meas = new Measurement();    
      meas->set_name(type);                         // the measurement type
      meas->set_hv(HyperVector(smear_m,meas_cov,type));  // the HyperVector 
      meas->set_name(RP::setup_volume,"mother");          // the volume
        
      // the true position of the measurement
      klog << log4cpp::Priority::DEBUG << "Create measurement plane " ;
      meas->set_position_hv(HyperVector(m,EMatrix(3,3,0),RP::xyz));

      //    Add the measurement to the vector
      fMeas.push_back(meas); 
    
      // ------ Store true information for pulls -----------
      if (betaTest == true) // we have a single electron of P=2.9 MeV
      {
        klog << log4cpp::Priority::DEBUG << "Storing true info " ;

        State state;
        TLorentzVector pos = beta->GetInitialVertex();
        TLorentzVector mom = beta->GetInitialMomentum();
        double q = beta->GetCharge();

        EVector v(fDim,0);
        for (int k=0;k<3;k++)
          v[k]=pos[k];
    
        v[3]=mom[0]/mom[2];
        v[4]=mom[1]/mom[2];

        klog << log4cpp::Priority::DEBUG << "model " << fModel << " dim " << fDim ;
        klog << log4cpp::Priority::DEBUG << " state vector " << v ;

        if (fKfs->Model() == "helix") 
          v[5]=q/(mom.Vect().Mag());
    
        state.set_hv(HyperVector(v,EMatrix(fDim,fDim,0),RP::slopes_curv_z));

        // create the node
        klog << log4cpp::Priority::DEBUG << " create the node "  ;
        Node* node = new Node();
        node->set_measurement(*meas);
        node->set_state(state);
        fSimTraj.nodes().push_back(node); 
      }
       
   }
    
    // add the measurements to the trajectory
    klog << log4cpp::Priority::INFO << " Measurement vector created "  ;
    Trajectory trj;
    trj.add_constituents(fMeas);  

    klog << log4cpp::Priority::INFO << " Trajectory --> " << trj ;

    return trj;
  }

//--------------------------------------------------------------------
  void KFSvcManager::SeedState(const ISvc& isvc, 
                        const Trajectory& traj, State& state) 
//--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    std::vector<const irene::Particle*> vBeta = isvc.Electrons();

    const irene::Particle* beta = vBeta[0];
    //auto betaTest = isvc.SelectBetaP(beta);

    EVector v(fDim,0);
    EMatrix C(fDim,fDim,0);

    klog << log4cpp::Priority::INFO << " ++KFSvcManager::SeedState --> " ;
    klog << log4cpp::Priority::INFO << " beta --> " 
          << isvc.PrintElectron(beta);

    const RecObject& meas = traj.nodes()[0]->measurement();
    const EVector& m = meas.vector();

    klog << log4cpp::Priority::INFO << " first measurement --> " 
          << m;

    TLorentzVector ul = beta->GetInitialMomentum();
    EVector u(3,0);
    u[0] = ul.X();
    u[1] = ul.Y();
    u[2] = ul.Z();

    klog << log4cpp::Priority::INFO << " Momentum vector --> " 
          << u;

    v[0] = m[0];   // a guess of the x position
    v[1] = m[1];   // a guess of the y position
    v[2] = m[2];  // a guess of the z position  
    v[3] = u[0]/u[2];  // a guess of dx/dz (use truth for the moment)
    v[4] = u[1]/u[2];  // a guess of dy/dz (use truth for the moment)

    klog << log4cpp::Priority::INFO << " Initial state --> " 
          << v;

    if (fKfs->Model() =="helix")
    {
      v[5]=-1./PMAX;  // assume forward fit!
    }
  else
    {
    // For the Straight line model q/p is a fix parameter use for MS computation 
    
    double qoverp=-1./PMAX;
    klog << log4cpp::Priority::INFO << " Sline model, set qoverp --> " << qoverp ;
    state.set_hv(RP::qoverp,HyperVector(qoverp,0));
    }

    // give a large diagonal covariance matrix
    C[0][0]= C[1][1]=100.*cm;
    C[2][2]= EGeo::zero_cov()/2.; // no error on Z since planes are perpendicular to z
    C[3][3]= C[4][4]=1;

    klog << log4cpp::Priority::INFO << " Cov matrix --> " << C ;

    if (fKfs->Model() =="helix")
      C[5][5]= 0.1;
  
    // Set the main HyperVector
    klog << log4cpp::Priority::INFO << " Set the main HyperVector --> "  ;
    state.set_hv(HyperVector(v,C));

    // Set the sense HyperVector (1=increasing z)
    double sense=u[2]/fabs(u[2]); // use truth for the moment

    klog << log4cpp::Priority::INFO << " Set the sense HyperVector --> " << sense ;
    state.set_hv(RP::sense,HyperVector(sense,0));

    // Set the model name
    klog << log4cpp::Priority::INFO << " Set the model name --> " << fModel ;
    state.set_name(fModel);

  // Set the representation (x,y,z dx/dz, dy/dz, q/p)
    klog << log4cpp::Priority::INFO << " Set the representation --> "  ;
    if (fKfs->Model() =="helix")
      state.set_name(RP::representation,RP::slopes_curv_z);
    else
      state.set_name(RP::representation,RP::slopes_z);
  }

//--------------------------------------------------------------------
  bool KFSvcManager::Fit(const ISvc& isvc) 
//--------------------------------------------------------------------
  { 
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::INFO << "In KFSvcManager::Fit --Creating Trajectory " ;
    auto traj = CreateTrajectory(isvc) ;

    klog << log4cpp::Priority::INFO << " Creating seed state "  ;
    State seed;
    SeedState(isvc,traj,seed); 
                        
      
      // fit the trajectory provided a seed state

    klog << log4cpp::Priority::INFO << " Fitting "  ;
      bool status = fRPMan.fitting_svc().fit(seed,traj);
      return status;

      // compute the pulls
      //      compute_pulls(_sim_traj,traj);

      //fill_histos(traj);
      
  }
//--------------------------------------------------------------------
  int KFSvcManager::ModelDim() const 
//--------------------------------------------------------------------
  {
    int dim;
    if (KFSetup::Instance().Model() =="helix")
    {
      dim = 6;                     // dimension of the state vector for this model
    } 
    else if (KFSetup::Instance().Model() =="sline")
    {
      dim = 5; 
    }
    else
    {
      std::cout << "ERROR: model unknown" << std::endl;
      exit (EXIT_FAILURE);
    }
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
}
