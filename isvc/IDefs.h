#ifndef IDEFS_
#define IDEFS_
/*
 typedefs for Irene
*/

#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <map>

#include <irene/Event.h>
#include <irene/Track.h>
#include <irene/Particle.h>
#include <TLorentzVector.h>

namespace alex {
typedef std::pair<TLorentzVector,double> IHit;
typedef std::vector<std::pair<TLorentzVector,double> > IHits;
typedef const irene::Particle* IParticle;
typedef std::vector<const irene::Particle*> IParticles;

}
#endif