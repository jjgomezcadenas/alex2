
// ----------------------------------------------------------------------------
//  $Id: ToyData.cpp 
//
//  Author:  <gomez@mail.cern.ch>
//  Created: July 2014
//  
//  Copyright (c) 2014 NEXT Collaboration
// ---------------------------------------------------------------------------- 



#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <map>
#include <iostream>
#include <sstream>

#include <alex/DParam.h>


using std::string; 
using std::cout;
using std::cin; 
using std::endl; 
using std::ostream;
using std::ostringstream;
using std::vector;

namespace alex {


//--------------------------------------------------------------------
  std::string DParam::Serialize() const
//--------------------------------------------------------------------
  {
    ostringstream s;
    
    s << "{Param: name = " << Name() << " Type = " << DataType()
        << " Value = " << Value() <<"}"<<endl <<std::ends;
    return s.str();
  }

}

