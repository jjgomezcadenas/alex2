
// ----------------------------------------------------------------------------
//  $Id: AXml.cc 
//
//  Author:  <gomez@mail.cern.ch>
//  Created: July 2014
//  
//  Copyright (c) 2014 NEXT Collaboration
// ---------------------------------------------------------------------------- 

#include <alex/AlexConfigure.h>
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
using std::pair;


using namespace tinyxml2;
namespace alex {

//--------------------------------------------------------------------
  void AlexConf::Init(std::string debugLevel, std::string rootName)  
//--------------------------------------------------------------------
  {
    fDebugLevel = debugLevel;
    fRootName = rootName;

    SetDebugLevel(fDebugLevel);
  }
//--------------------------------------------------------------------
  void AlexConf::ParseConfiguration(std::string configFile)
//--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::INFO << " config file =" << configFile;

    fDoc.LoadFile( configFile.c_str() );
    if (fDoc.ErrorID()!=0) 
    {
      klog << log4cpp::Priority::ERROR << "Failed loading config file error = " << fDoc.ErrorID();
      exit (EXIT_FAILURE);
    }

    std::pair<string,string> tags;
    tags.first = "path";
    tags.second="name";

    XMLElement* rootElement =ParseRoot(); 

    const XMLElement* algoElement = rootElement->FirstChildElement ("Algos") ;
    klog << log4cpp::Priority::DEBUG << " algoElement " << algoElement->Name();
    std::pair<string,string> algoPathName =  ParseStringPair(algoElement,tags);

    const XMLElement* dstElement = algoElement->NextSiblingElement ("DST") ;
    klog << log4cpp::Priority::DEBUG << " dstElement " << dstElement->Name();
    std::pair<string,string> dstPathName =  ParseStringPair(dstElement,tags);

    const XMLElement* histoElement = algoElement->NextSiblingElement ("HistoFile") ;
    klog << log4cpp::Priority::DEBUG << " histoElement " << histoElement->Name();
    std::pair<string,string> histoPathName =  ParseStringPair(histoElement,tags);

    const XMLElement* eventElement = algoElement->NextSiblingElement ("Events") ;
    klog << log4cpp::Priority::DEBUG << " eventElement " << eventElement->Name();

    tags.first = "runMax";
    tags.second="runDebug";

    std::pair<int,int> events =  ParseIntPair(eventElement,tags);

    const XMLElement* debugElement = algoElement->NextSiblingElement ("Debug") ;
    klog << log4cpp::Priority::DEBUG << " debugElement " << debugElement->Name();

    const XMLElement* elem = debugElement->FirstChildElement ("level") ;
    string debug  = elem->GetText();
    klog << log4cpp::Priority::DEBUG << " debug = " << debug;

  }

//--------------------------------------------------------------------
  XMLElement* AlexConf::ParseRoot() 
//--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();

    XMLElement* rootElement = fDoc.RootElement();
    string rootName = rootElement->Name();
    if (rootName != fRootName)
    {
      klog << log4cpp::Priority::ERROR << "Configuration file root should be <" <<fRootName <<">";
      exit (EXIT_FAILURE);
    }
    return rootElement;
  }

//--------------------------------------------------------------------
  std::pair<string,string> AlexConf::ParseStringPair(const XMLElement* mom, 
                                                    const std::pair<string,string>& tags) 
//--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();

    const XMLElement* elem = mom->FirstChildElement (tags.first.c_str()) ;
    string first = elem->GetText();
    klog << log4cpp::Priority::DEBUG << " first = " << first;

    const XMLElement*  nextElem = elem->NextSiblingElement (tags.second.c_str()) ;
    string second = nextElem->GetText();
    klog << log4cpp::Priority::DEBUG << " second = " << second;

    std::pair<string,string> pathName;
    pathName.first = first;
    pathName.second = second;
    return pathName;

  }

//--------------------------------------------------------------------
  std::pair<int,int> AlexConf::ParseIntPair(const XMLElement* mom, 
                                               const std::pair<string,string>& tags) 
//--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();

    const XMLElement* elem = mom->FirstChildElement (tags.first.c_str()) ;
    string first = elem->GetText();
    klog << log4cpp::Priority::DEBUG << " first = " << first;

    const XMLElement*  nextElem = elem->NextSiblingElement (tags.second.c_str()) ;
    string second = nextElem->GetText();
    klog << log4cpp::Priority::DEBUG << " second = " << second;

    std::pair<int,int> intPair;
    intPair.first = atoi(first.c_str());
    intPair.second = atoi(second.c_str());
    return intPair;

  }


}

