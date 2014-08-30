
/*
 Hit
 JJGC, July, 2014.
*/

#include "Hit.h"


namespace alex {


	Hit::Hit(const IHit&  ihit)
	{
		auto xyz =ihit.first;
		fEne =ihit.second;
		fXYZ(0)=xyz.X();
		fXYZ(1)=xyz.Y();
		fXYZ(2)=xyz.Z();
	}
	Hit::Hit(const Hit& hit)
	{
		fEne = hit.Edep();
		fXYZ = hit.XYZ();
		fType = hit.Type();
		fOwner = hit.Owner();
	}
		

}
