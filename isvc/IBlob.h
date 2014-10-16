#ifndef IBLOB_
#define IBLOB_
/*
 IBLOB represents a blob at the end of the IBeta track.
 A blob is a set of hits 
 JJGC, JR, October, 2014.

 Notice: the sorted hit vector should sort the hits of the blob in z
 this implies that the first "true hit" from the electron track could be not the first
*/

#include "Hit.h"

namespace alex {

class IBlob {
      public:
              IBlob(){};
              IBlob(const IBlob&  iblob);
              virtual ~IBlob();
              void AddHit(const Hit& hit);
              void AddHit(const IHit& hit);
              double GetEnergy() const;
       
              void CreateBlob(double x, double y, double z, double E);
               
              std::vector<const Hit* > GetHit() const {return fHitVector;}
              std::vector<const Hit* > GetSortedHits() const {return fSortedHitVector;}
              const Hit* GetBlob() const {return fBlob;}

      private:
              void Clear();
              void DeleteHits();
              std::vector<const Hit*> fHitVector;
              std::vector<const Hit*> fSortedHitVector;  //hits sorted (normally in z)
              Hit* fBlob;
                      
      };

}
#endif
