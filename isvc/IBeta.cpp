
/*
 Hit
 JJGC, July, 2014.
*/

#include "IBeta.h"


namespace alex {

	IBeta::IBeta(const IBeta& ibeta)
	{
		fType = ibeta.GetType();
		fMcType =ibeta.GetMCType();
		fV0 = ibeta.GetV0(); 
		fVx =ibeta.GetVX();
	}
	void IBeta::SetVX(TVector3 V1, TVector3 V2) 
		{
			fVx.first = V1;
			fVx.second = V2;
		}	

}
