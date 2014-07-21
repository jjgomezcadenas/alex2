
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

#include <alex/LogUtil.h>
#include <alex/ToyData.h>

using std::string; 
using std::cout;
using std::cin; 
using std::endl; 
using std::ostream;
using std::ostringstream;
using std::vector;

namespace alex {


//--------------------------------------------------------------------
  std::string ToyData::Serialize() const
//--------------------------------------------------------------------
  {
    ostringstream s;
   
    //s << std::endl;
    s<<"<DataBucket>";
    s << "\t<name>"<<this->Name()<<"</name>"<<endl;
    s<<"\t<type>TVector3</type>"<<endl;
    s<<"\t<dim>"<<3<<"</dim>"<<endl;
    s<<"\t<value>"<<fX[0] <<" " <<fX[1] <<" " <<fX[2] <<"</value>" <<endl; 
    s<<"</DataBucket>" <<endl <<std::ends;
    return s.str();

  }
//--------------------------------------------------------------------
  void ToyData::Recreate(std::string)
//--------------------------------------------------------------------
  {
  }


}

