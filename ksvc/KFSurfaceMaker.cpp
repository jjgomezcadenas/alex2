#include "KFSurfaceMaker.h"
#include <recpack/HelixEquation.h>
#include <recpack/Plane.h>

namespace alex {

//*******************************************
  KFSurfaceMaker::KFSurfaceMaker(MatchingSvc* mat, Setup* setup) : ISurfaceMaker(){
//*******************************************

    std::cout << "Set up KFSurfaceMaker." << std::endl;
    _matching = mat;
    _setup = setup;

  }


//*******************************************
  bool KFSurfaceMaker::get_surface(const Measurement& meas, const State& state, Surface& surf) {
//*******************************************

    // Get the position of the measurement.
    const HyperVector& pos_hv = meas.position_hv();

    // Get the direction of the state.
    EVector u = _equation->direction(state);

    if (verbosity(VVERBOSE)) {
      std::cout <<  "  KFSurfaceMaker::get_surface(). (for meas) u = " << print(u) << ", pos = " << print(pos_hv.vector()) << std::endl;    
    }    
 
    return get_surface(meas, pos_hv, state, surf);
  }

//*******************************************
  bool KFSurfaceMaker::get_surface(const Node& node, const State& state, Surface& surf) {
//*******************************************
  
    const Measurement& meas = node.measurement();

    // Get the direction of the state.
    EVector u = _equation->direction(state);

    // Get the position of the node.
    EVector pos = _equation->position(node.state());
    HyperVector pos_hv(pos,EMatrix(3,3,0));

    return get_surface(meas, pos_hv, state, surf);
  }

//*******************************************
  bool KFSurfaceMaker::get_surface(const Measurement& meas, const HyperVector& pos_hv, const State& state, Surface& surf) {
//*******************************************
  
    EVector u = _equation->direction(state);

    /*EVector axis(3,0);

    // Set the surface normal according to the z-direction.
    if(u[2] < 0) { 
      axis[2] = -1.;
      std::cout << "Negative normal." << std::endl;
    }
    else {
      axis[2] = 1.;
      std::cout << "Positive normal." << std::endl;
    }

    // Create the surface.
    surf = Plane(pos_hv,axis);*/

    EVector axis(3,0);

    // Set the axis appropriately if we have failed in propagating to any of the surfaces at this point.
    if(state.name("failextrap") == "xy") axis[2] = 1;
    else if(state.name("failextrap") == "xz") axis[1] = 1;
    else if(state.name("failextrap") == "yz") axis[0] = 1;
    else {

      // Ensure we don't choose a normal in a direction in which the propagation already failed.
      int failAxis = -2;
      if(state.name("failextrap") == "x") failAxis = 0;
      else if(state.name("failextrap") == "y") failAxis = 1;
      else if(state.name("failextrap") == "z") failAxis = 2;

      // Find the directional component with the largest value.
      double u_max = 0;
      int i_max = -1;
      for(int i = 0; i < 3; i++) {
        if ((fabs(u[i]) >= u_max || i_max < 0) && i != failAxis) {
          u_max = fabs(u[i]);
          i_max = i;
        }    
      }

      if (verbosity(VVERBOSE)) {
        std::cout <<  "    --> i_max = " << i_max << ", u[i_max], u_max = " << u[i_max] << ", " << u_max << std::endl;
      }

      // Make sure that measurement has an associated volume.
      if (!meas.names().has_key(RP::setup_volume)) {
        if (verbosity(VERBOSE))
          std::cout<<"KFSurfaceMaker::get_surface(): not volume for measurement "<<meas<<std::endl;
        return false;
      }

      // Choose the largest directional component as the axis.
      axis[i_max] = 1;
    }

    // Create the surface.
    surf = Plane(pos_hv,axis);

    return true;  
  }

}
