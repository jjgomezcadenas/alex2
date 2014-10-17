
/*
 Hit
 JJGC, July, 2014.
*/

#include "Hit.h"


namespace alex {


	Hit::Hit(const IHit&  ihit)
	{
		auto xyz =ihit.first;
		fEne = ihit.second;
		fXYZ(0) = xyz.X();
		fXYZ(1) = xyz.Y();
		fXYZ(2) = xyz.Z();
                fMinDist = -1.;
	}
	Hit::Hit(const Hit& hit)
	{
		fEne = hit.Edep();
		fXYZ = hit.XYZ();
		fType = hit.Type();
		fOwner = hit.Owner();
                fMinDist = hit.MinDist();
	}
	
	Hit::Hit(double x, double y, double z, double E)
	{
		fEne = E;
		fXYZ = TVector3(x,y,z);
		fType = "BLOB";
		fOwner = 0;
                fMinDist = -1.;
	}

	Hit::Hit(double x, double y, double z, double zend, double E)
	{
		fEne = E;
		fXYZ = TVector3(x,y,z);
		fType = "EFF";
		fOwner = 0;
		fMinDist = -1;
	}
}
