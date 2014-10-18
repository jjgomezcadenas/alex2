#ifndef KFRP_
#define KFRP_
/*
 Kalman Filter RecPack
 

 JJGC, October, 2014.
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


/*
RECPACK

Vector and Matrices: EVector and EMatrix are just typedefs: RP uses CLHEP vector and matrices

class HyperVector
------------------

An HyperVector is a agregation of a vector and a covariance-matrix
found in recpack/util

protected:
    
    //! vector
    EVector _vector;
    
    //! covariance matrix
    EMatrix _matrix;

    //! the hypervector representation
    std::string _rep; 

    //! degrees of freedom
    int _ndof;

    //! the chi2 of the HyperVector
    double _chi2;

An HyperVector Object
 - Has a main HyperVector (a vector and a covariance matrix)
    - Optional: Has a dictionary of extra HyperVectors

A NamedObject has a name and a dictionary of names

protected:

    //! name of the NamedObject
    Key _name;
      
    //! dictionary of names
    dict::dictionary<Key> _names;

------
A State
--------
  
    - State is a NamedObject
    - State is a HyperVectorObject
    - State has a running parameter  

-------
A RecResult
-------
  
    - RecResult has a dictionary of status (bool)
    - RecResult has a quality (double), i.e chi squared
    - RecResult has a number of degrees of freedom
    - Optional: RecResult has a dictionary of extra qualities (doubles)

-----
 Node
-----
   
    A Node is the agregation of a state and a RecObject. A collection of 
    residuals and the quality of agreement between the state projection and
    the RecObject.
    - Node is a NamedObject
    - Node is a RecResult
    - Node has a local state
    - Node has a RecObject
    - Node has a residual HyperVectorObject
    _ Node has a surface in which the residuals are computed
 
     //! local state 
    State* _state;

    //! RecObject 
 
    RecObject* _RecObject;

    //! hypervectors for  residuals
    HyperVectorObject _residuals;

    //! surface where the residual is computed
    const Surface* _surface;

    

-------
RecObject
---------
  
    A RecObject is a collection of Nodes: a collection of measurements and
    states with an agreement between them.
    - RecObject is a RecResult 
    - RecObject is a NamedObject
    - RecObject is a HyperVectorObject
    - RecObject has a collection of nodes.
    - RecObject has a collection of measurements

------
A Trajectory
------
    A Trajectory is a collection of Nodes: a collection of RecObjects (measurements) and
    states with an agreement between them.
    - Trajectory is a RecObject 
    - Trajectory has a length
 
  class Trajectory: public RecObject


----
A Measurements and a measurement vector: They are RecObjects
typedef RecObject Measurement;

typedef std::vector<Measurement*> measurement_vector;



-----------
Regarding which info can  be retrieved from nodes as a result of the
kalman filtering:

different state hypervector values and residuals are stored with the
corresponding names:

E.g. to retrieve predicted dynamic vector and matrix from prediction step:
// get the state
State& state = node.state();
const HyperVector& hv_P = state.hv(RP::predicted);
const EVector& x_P = hv_P.vector();
const EMatrix& C_P = hv_P.matrix();

// Retrieve predicted residual
const HyperVector resHV_P = node.residuals().hv(predicted_name);
const EVector& r_P = resHV_P.vector();
const EMatrix& R_P = resHV_P.matrix();

node.quality(predicted_name,chi2); //or resHV_P.chi2()

available names:
RP:predicted, RP::filtered
RP::smoothed (default)

RP::predicted_back and RP::weighted when using back filtering isntead of
smoothing

add_projector(RP::pos, new ComponentProjector(RP::pos_dir_curv, 0,1,2));
    add_projector(RP::dir, new ComponentProjector(RP::pos_dir_curv, 3,4,5));
    add_projector(RP::curv, new ComponentProjector(RP::pos_dir_curv, 6) );
    add_projector(RP::pos_dir, new ComponentProjector(RP::pos_dir_curv, 0,1,2,3,4,5) );
    add_projector(RP::pos_dir_curv, new ComponentProjector(RP::pos_dir_curv, 0,1,2,3,4,5,6) );

    add_projector(RP::slope, new ComponentProjector(RP::slopes_curv, 3,4) );
    add_projector(RP::pos_slope, new ComponentProjector(RP::slopes_curv, 0,1,2,3,4) );
    add_projector(RP::slopes_curv,   new ComponentProjector(RP::slopes_curv,   0,1,2,3,4,5) );
    add_projector(RP::slopes_curv_x, new ComponentProjector(RP::slopes_curv_x, 0,1,2,3,4,5) );
    add_projector(RP::slopes_curv_y, new ComponentProjector(RP::slopes_curv_y, 0,1,2,3,4,5) );
    add_projector(RP::slopes_curv_z, new ComponentProjector(RP::slopes_curv_z, 0,1,2,3,4,5) );

    add_projector(RP::xyz, new ComponentProjector(RP::pos_dir_curv, 0,1,2) );
    add_projector(RP::xy, new ComponentProjector(RP::pos_dir_curv, 0,1) );
    add_projector(RP::xz, new ComponentProjector(RP::pos_dir_curv, 0,2) );
    add_projector(RP::yz, new ComponentProjector(RP::pos_dir_curv, 1,2) );
    add_projector(RP::x, new ComponentProjector(RP::pos_dir_curv, 0) );
    add_projector(RP::y, new ComponentProjector(RP::pos_dir_curv, 1) );
    add_projector(RP::z, new ComponentProjector(RP::pos_dir_curv, 2) );
    add_projector(RP::time, new ComponentProjector(RP::pos_dir_curv_time, 7) );
  
*/


#endif