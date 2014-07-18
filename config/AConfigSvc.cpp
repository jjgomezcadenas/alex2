
// ----------------------------------------------------------------------------
//  $Id: AXml.cc 
//
//  Author:  <gomez@mail.cern.ch>
//  Created: July 2014
//  
//  Copyright (c) 2014 NEXT Collaboration
// ---------------------------------------------------------------------------- 

#include <alex/AConfigSvc.h>
#include <alex/LogUtil.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>
#include <sstream>
#include <cstdio>
#include <cstring>

using std::string; 
using std::cout; 
using std::endl; 
using std::ostream;
using std::vector;


using namespace tinyxml2;
namespace alex {


//--------------------------------------------------------------------
  void AXml::ParseParamFile(std::string paramFile)
//--------------------------------------------------------------------
  {
    
    SetDebugLevel("DEBUG");
    log4cpp::Category& klog = log4cpp::Category::getRoot();


    fDoc.LoadFile( paramFile.c_str() );
    if (fDoc.ErrorID()!=0) 
    {
      klog << log4cpp::Priority::ERROR << "Failed loading param file";
      exit (EXIT_FAILURE);
    }
    
    XMLElement* rootElement = fDoc.RootElement();
    klog << log4cpp::Priority::DEBUG << "name of root element" << rootElement->Name();


    const XMLElement* param = rootElement->FirstChildElement ("Param") ;
    klog << log4cpp::Priority::DEBUG << "FirstChildElement (must be Param) " << param->Name();

    ParseParamElement(param);

    const XMLElement*  nextParam = param->NextSiblingElement ("Param") ;

    while (nextParam != NULL)
    {
      klog << log4cpp::Priority::DEBUG << "Next Sibling Param (must be Param) " << nextParam->Name();
      ParseParamElement(nextParam);
      nextParam = nextParam->NextSiblingElement ("Param") ;
    }

    const XMLElement* array = rootElement->FirstChildElement ("Array") ;
    ParseArrayElement(array);

    const XMLElement*  nextArray = array->NextSiblingElement ("Array") ;

    while (nextArray != NULL)
    {
      ParseParamElement(nextArray);
      nextArray = nextArray->NextSiblingElement ("Array") ;
    }
    

  }
//--------------------------------------------------------------------
  void AXml::ParseParamElement(const XMLElement* param) const
//--------------------------------------------------------------------
  {
    //   <Param>
    //   <name>dataPath</name>
    //   <dataType>string</dataType>
    //   <value>/Users/jjgomezcadenas/Development/NEXT/DATA</value>
    // </Param>
    
    log4cpp::Category& klog = log4cpp::Category::getRoot();

    const XMLElement* nameParam = param->FirstChildElement ("name") ;
    string textNameParam = nameParam->GetText();
    klog << log4cpp::Priority::DEBUG << " Param name text = " << textNameParam;

    const XMLElement* dataTypeParam = param->FirstChildElement ("dataType") ;
    string dataTypeParamText= dataTypeParam->GetText();
    klog << log4cpp::Priority::DEBUG << "Param data type text =" << dataTypeParamText;

    const XMLElement* ValueParam = param->FirstChildElement ("value") ;
    string ValueParamText= ValueParam->GetText();
    klog << log4cpp::Priority::DEBUG << "Param value text =" << ValueParamText;
  }
//--------------------------------------------------------------------
  void AXml::ParseArrayElement(const XMLElement* array) const
//--------------------------------------------------------------------
  {
    //   <Array>
    //   <name>P</name>
    //   <dataType>double</dataType>
    //   <dim>3</dim>
    //   <value>0.1 1.0 3</value>
    // </Array>

    log4cpp::Category& klog = log4cpp::Category::getRoot();

    const XMLElement* name = array->FirstChildElement ("name") ;
    string textName = name->GetText();
    klog << log4cpp::Priority::DEBUG << " Array name text = " << textName;

    const XMLElement* dataType = array->FirstChildElement ("dataType") ;
    string textDataType= dataType->GetText();
    klog << log4cpp::Priority::DEBUG << "Array data type text =" << textDataType;

    const XMLElement* dim = array->FirstChildElement ("dim") ;
    string textDim= dim->GetText();
    klog << log4cpp::Priority::DEBUG << "Array dim text =" << textDim;

    const XMLElement* value = array->FirstChildElement ("value") ;
    string textValue= value->GetText();
    klog << log4cpp::Priority::DEBUG << "Array value text =" << textValue;
  }

}

