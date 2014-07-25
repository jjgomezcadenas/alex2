
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
#include <tinyxml2.h>

using std::string; 
using std::cout;
using std::cin; 
using std::endl; 
using std::ostream;
using std::ostringstream;
using std::vector;

using namespace tinyxml2;
namespace alex {


//--------------------------------------------------------------------
  std::string ToyData::Serialize() const
//--------------------------------------------------------------------
  {
    ostringstream s;
   
    //s << std::endl;
    s<<"<DataBucket>\n";
    s << "\t<name>"<<this->Name()<<"</name>"<<endl;
    s<<"\t<type>TVector3</type>"<<endl;
    s<<"\t<dim>"<<3<<"</dim>"<<endl;
    s<<"\t<value>"<<fX[0] <<" " <<fX[1] <<" " <<fX[2] <<"</value>" <<endl; 
    s<<"</DataBucket>" <<endl <<std::ends;
    return s.str();

  }
// //--------------------------------------------------------------------
//   void ToyData::Recreate(std::string xml)
// //--------------------------------------------------------------------
//   {
//     log4cpp::Category& klog = log4cpp::Category::getRoot();

//     XMLDocument doc;
//     doc.Parse( xml.c_str() );
//     if (doc.ErrorID()!=0) 
//     {
//       klog << log4cpp::Priority::ERROR << " In ToyData::Recreate-- Failed parsing";
//       exit (EXIT_FAILURE);
//     }

//     XMLElement* rootElement = doc.RootElement();
//     klog << log4cpp::Priority::DEBUG << "name of root element " << rootElement->Name();


//     const XMLElement* firstElement = rootElement->FirstChildElement ("name") ;
//     klog << log4cpp::Priority::DEBUG << "First Element = " << firstElement->Name();
//     string text = firstElement->GetText();
//     klog << log4cpp::Priority::DEBUG << " text = " << text;

//     const XMLElement* nextElement = firstElement->NextSiblingElement ("type") ;
//     klog << log4cpp::Priority::DEBUG << " Next Element = " << nextElement->Name();
//     text = nextElement->GetText();
//     klog << log4cpp::Priority::DEBUG << " text = " << text;

//     nextElement = firstElement->NextSiblingElement ("dim") ;
//     klog << log4cpp::Priority::DEBUG << " Next Element = " << nextElement->Name();
//     text = nextElement->GetText();
//     klog << log4cpp::Priority::DEBUG << " text = " << text;

//     nextElement = firstElement->NextSiblingElement ("value") ;
//     klog << log4cpp::Priority::DEBUG << " Next Element = " << nextElement->Name();
//     text = nextElement->GetText();
//     klog << log4cpp::Priority::DEBUG << " text = " << text;

//   }

}

