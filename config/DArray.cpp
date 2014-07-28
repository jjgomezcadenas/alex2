
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

#include <alex/DArray.h>
#include <alex/StringOperations.h>


using std::string; 
using std::cout;
using std::cin; 
using std::endl; 
using std::ostream;
using std::ostringstream;
using std::vector;

namespace alex {


//--------------------------------------------------------------------
  std::string DArray::Serialize() const
//--------------------------------------------------------------------
  {
    ostringstream s;
    
    s << "{Array: name = " << Name() << " Type = " << DataType()
        << " Dim = "<< Dim() <<  " Value = " << Value() <<"}"<<endl <<std::ends;
    return s.str();
  }

//--------------------------------------------------------------------
  std::vector<std::string>  DArray::Split() 
//--------------------------------------------------------------------
  {
  	return SplitString(Value());
  }

}

