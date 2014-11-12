/* -*- mode: c++ -*- */
#ifndef KFSurfaceMaker_
#define KFSurfaceMaker_
/*
 *  KFSurfaceMaker: creates 
 *   
 *
 *    Josh, October, 2014.
 *    Based on TRecPackMatchingSurfaceMaker
 *    */

#include <alex/KFRecpack.h>

namespace alex {

class KFSurfaceMaker : public ISurfaceMaker {

  public:
    
    KFSurfaceMaker(MatchingSvc* mat, Setup* setup);

    //! default destructor
    virtual ~KFSurfaceMaker() {};

    //!  Gets the surface provided the measurement
    virtual bool get_surface(const Measurement& meas, const State& state, Surface& surf);

    //!  Gets the surface provided the node
    virtual bool get_surface(const Node& node, const State& state, Surface& surf);

    //!  Basic function to get the surface
    virtual bool get_surface(const Measurement& meas, const HyperVector& pos_hv, const State& state, Surface& surf);

  protected:

    MatchingSvc* _matching;
    Setup* _setup;

  };

} 
#endif
