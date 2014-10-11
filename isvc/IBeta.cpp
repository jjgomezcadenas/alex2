
/*
 Hit
 JJGC, July, 2014.
*/

#include "IBeta.h"


namespace alex {

	void IBeta::AddBlob(const IBlob&  blob)
	{
		auto iblob =new IBlob(blob);
		fBlobs.push_back(iblob);
	}

  void IBeta::ClearBlobs()
  {
    fBlobs.clear();
  }

	IBeta::~IBeta()
	{
		for(auto iblob : fBlobs)
			delete iblob; 

		ClearBlobs();	
	}
}
