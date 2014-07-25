
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

#include <alex/DHistos.h>


using std::string; 
using std::cout;
using std::cin; 
using std::endl; 
using std::ostream;
using std::ostringstream;
using std::vector;

namespace alex {


//--------------------------------------------------------------------
  std::string DH1::Serialize() const
//--------------------------------------------------------------------
  {
    ostringstream s;
    
    s << "{H1D: name = " << Name() << " Title = " << Title()
        << " nbinsx = "<< Nbinsx() 
        <<  " Xlow = " << Xlow() <<  " Xup = " << Xup() <<"}"<< endl <<std::ends;
    return s.str();
  }

//--------------------------------------------------------------------
  std::string DH2::Serialize() const
//--------------------------------------------------------------------
  {
    ostringstream s;
    
    s << "{H2D: name = " << Name() << " Title = " << Title()
        << " nbinsx = "<< Nbinsx() 
        << " nbinsy = "<< Nbinsy() 
        <<  " Xlow = " << Xlow() <<  " Xup = " << Xup()
        <<  " Ylow = " << Ylow() <<  " Yup = " << Yup() 
        <<"}"<< endl <<std::ends;
    return s.str();
  }
}

