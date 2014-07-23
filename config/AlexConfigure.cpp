
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
#include <alex/StringOperations.h>
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
using std::ostringstream;


using namespace tinyxml2;
namespace alex {

//--------------------------------------------------------------------
  void AlexConf::Init(std::string debugLevel, std::string rootName)  
//--------------------------------------------------------------------
  {
    fDebugLevel = debugLevel;
    fRootName = rootName;

    
    fStags.first = "path";
    fStags.second="name";

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
      klog << log4cpp::Priority::ERROR 
      << "ParseConfiguration::Failed loading config file error = " 
      << fDoc.ErrorID();
      exit (EXIT_FAILURE);
    }


    XMLElement* rootElement = fDoc.RootElement();
    klog << log4cpp::Priority::DEBUG << " rootElement " << rootElement->Name();

    const XMLElement* algoElement = rootElement->FirstChildElement ("Algos") ;
    klog << log4cpp::Priority::DEBUG << " algoElement " << algoElement->Name();
    fAlgosPathName =  ParseStringPair(algoElement,fStags);

    const XMLElement* dstElement = algoElement->NextSiblingElement ("DST") ;
    klog << log4cpp::Priority::DEBUG << " dstElement " << dstElement->Name();
    fDstPathName =  ParseStringPair(dstElement,fStags);

    const XMLElement* histoElement = algoElement->NextSiblingElement ("HistoFile") ;
    klog << log4cpp::Priority::DEBUG << " histoElement " << histoElement->Name();
    fHistoPathName =  ParseStringPair(histoElement,fStags);

    const XMLElement* eventElement = algoElement->NextSiblingElement ("Events") ;
    klog << log4cpp::Priority::DEBUG << " eventElement " << eventElement->Name();

    std::pair<std::string,std::string> tags;
    tags.first = "runMax";
    tags.second="runDebug";

    fEvents =  ParseIntPair(eventElement,tags);

    const XMLElement* debugElement = algoElement->NextSiblingElement ("Debug") ;
    klog << log4cpp::Priority::DEBUG << " debugElement " << debugElement->Name();

    const XMLElement* elem = debugElement->FirstChildElement ("level") ;
    string fDebug  = elem->GetText();
    klog << log4cpp::Priority::DEBUG << " debug = " << fDebug;

    ParseAlgosConfiguration();

    ParseAlgos();

  }
//--------------------------------------------------------------------
  void AlexConf::ParseAlgos()
//--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::DEBUG << " In ParseAlgos() " ;

    for (auto algoPath : fAlgoPath)
    {
      klog << log4cpp::Priority::DEBUG << " algoPath = " << algoPath;
      fDoc.LoadFile( algoPath.c_str() );
      if (fDoc.ErrorID()!=0) 
      {
        klog << log4cpp::Priority::ERROR 
        << "ParseConfiguration::Failed loading config file error = " 
        << fDoc.ErrorID();
        exit (EXIT_FAILURE);
      }

      XMLElement*  rootElement = fDoc.RootElement();
      klog << log4cpp::Priority::DEBUG << " rootElement " << rootElement->Name();

      const XMLElement* param = rootElement->FirstChildElement ("Param") ;
      klog << log4cpp::Priority::DEBUG << "FirstChildElement (must be Param) " << param->Name();

      if (param != NULL)
      {
        ParseParamElement(param);
        param = param->NextSiblingElement ("Param") ;
        while (param != NULL)
        {
          klog << log4cpp::Priority::DEBUG << "Next Sibling Param (must be Param) " << param->Name();
          ParseParamElement(param);
          param = param->NextSiblingElement ("Param") ;
        }
      }

      const XMLElement* array = rootElement->FirstChildElement ("Array") ;

      if (array != NULL)
      {
        ParseArrayElement(array);
        array = array->NextSiblingElement ("Array") ;
        while (array != NULL)
        {
          ParseParamElement(array);
          array = array->NextSiblingElement ("Array") ;
        }
      }
    

//   <H1D>
//     <name>fHY</name>
//     <title>Y distribution</title>
//     <nbinsx>50</nbinsx>
//     <xlow>-10.0</xlow>
//     <xup>10.0</xup>
//   </H1D>
//   <H2D>
//     <name>fHYZ</name>
//     <title>Y vs Z distribution</title>
//     <nbinsx>10</nbinsx>
//     <nbinsy>10</nbinsy>
//     <xlow>-10.0</xlow>
//     <xup>10.0</xup>
//     <ylow>-10.</ylow>
//     <yup>10.0</yup>
//   </H2D>
// </ToyAnalysis2>
    }
  }
//--------------------------------------------------------------------
  void AlexConf::ParseAlgosConfiguration()
//--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    string xmlPath = PathFromStrings(fAlgosPathName.first,fAlgosPathName.second);

    klog << log4cpp::Priority::DEBUG << "ParseAlgosConfiguration:: xmlPath ="
         << xmlPath;
    
    //tinyxml2::XMLDocument doc;
    fDoc.LoadFile( xmlPath.c_str() );
    if (fDoc.ErrorID()!=0) 
    {
      klog << log4cpp::Priority::ERROR 
      << "ParseAlgosConfiguration::Failed loading config file error = " 
      << fDoc.ErrorID();
      exit (EXIT_FAILURE);
    }

    XMLElement* rootElement = fDoc.RootElement();
    klog << log4cpp::Priority::DEBUG << " rootElement " << rootElement->Name();

    const XMLElement* algoElement = rootElement->FirstChildElement ("Algo") ;
    klog << log4cpp::Priority::DEBUG 
    << " ParseAlgosConfiguration:: algoElement " << algoElement->Name();

    std::pair<std::string,std::string> algoPathName = ParseStringPair(algoElement,fStags);
    string algoPath = PathFromStrings(algoPathName.first,algoPathName.second);
    klog << log4cpp::Priority::DEBUG 
    << " algoPath " << algoPath;

    fAlgoNames.push_back(algoPathName.second);
    algoPath = MergeStrings(algoPath,".xml");
    fAlgoPath.push_back(algoPath);

    const XMLElement*  nextAlgo = algoElement->NextSiblingElement ("Algo") ;

    while (nextAlgo != NULL)
    {
      algoPathName = ParseStringPair(nextAlgo,fStags);
      algoPath = PathFromStrings(algoPathName.first,algoPathName.second);
      algoPath = MergeStrings(algoPath,".xml");
      klog << log4cpp::Priority::DEBUG 
      << " algoPath " << algoPath;

      fAlgoNames.push_back(algoPathName.second);
      fAlgoPath.push_back(algoPath);
      nextAlgo = nextAlgo->NextSiblingElement ("Algo") ;
    }

  }

//--------------------------------------------------------------------
  std::string AlexConf::SerializeAlgoNames() const
//--------------------------------------------------------------------
  {
    return SerializeVectorInList(fAlgoNames);
  }
//--------------------------------------------------------------------
  std::string AlexConf::SerializeAlgoPaths() const
//--------------------------------------------------------------------
  {
    return SerializeVectorInList(fAlgoPath);
  }
//--------------------------------------------------------------------
  std::string AlexConf::SerializeAConfHeader() const
//--------------------------------------------------------------------
  {
    ostringstream s;

    s<<"\n#ifndef AACONF_" << endl;
    s<<"#define AACONF_" << endl;
    s<<"// Generated by AlexConf: do not edit" << endl;

    s<<"#include <string>" << endl;
    s<<"#include <vector>" << endl;
    s<<"#include <utility>" << endl;
    s<<"#include <memory>" << endl;
    s<<"#include <map>" << endl;

    s<<"namespace alex {" << endl;

    s<<"  class AAConf {" << endl;
    s<<"    public:" << endl;
    s<<"      AAConf();" << endl;
    s<<"      virtual ~AAConf(){};" << endl;
    s<<"      std::string AlgosPath() const {return fAlgoPathName.first;}" << endl;
    s<<"      std::string AlgosName() const {return fAlgoPathName.second;}" << endl;
    s<<"      std::string DstPath() const {return fDstPathName.first;}" << endl;
    s<<"      std::string DstName() const {return fDstPathName.second;}" << endl;
    s<<"      std::string HistoPath() const {return fHistoPathName.first;}" << endl;
    s<<"      std::string HistoName() const {return fHistoPathName.second;}" << endl;
    s<<"      int EventsToRun()const {return fEvents.first;}" << endl;
    s<<"      int EventsToDebug()const {return fEvents.second;}" << endl;
    s<<"      std::string DebugLevel()const {return fDebug;}" << endl;
    s<<"    private:" << endl;
    
    s<<"      std::pair<std::string,std::string> fAlgoPathName;" << endl;
    s<<"      std::pair<std::string,std::string> fDstPathName;" << endl;
    s<<"      std::pair<std::string,std::string> fHistoPathName;" << endl;
    s<<"      std::pair<int,int> fEvents;" << endl;
    s<<"      std::string fDebug;" << endl;
    s<<"  };" << endl;
    s<<"}" << endl;
    s<<"#endif" << endl;
    return s.str();
  }

//--------------------------------------------------------------------
  std::string AlexConf::SerializeAConfCPP() const
//--------------------------------------------------------------------
  {
    ostringstream s;
    s<<"\n// Generated by AlexConf: do not edit" << endl;

    s<<"#include <AAConf.h>" << endl;

    s<<"AAConf::AAConf()" << endl;
    s<<"{" << endl;
    s<<"  fAlgoPathName.first=" << fAlgosPathName.first << endl;
    s<<"  fAlgoPathName.second=" << fAlgosPathName.second << endl;
    s<<"  fDstPathName.first=" << fDstPathName.first << endl;
    s<<"  fDstPathName.second=" << fDstPathName.second << endl;
    s<<"  fHistoPathName.first=" << fDstPathName.first << endl;
    s<<"  fHistoPathName.second=" << fDstPathName.second << endl;
    s<<"  fEvents.first=" << fEvents.first << endl;
    s<<"  fEvents.second=" << fEvents.second << endl;
    s<<"  fDebug=" << fDebug << endl;
    s<<"}" << endl;
    
    return s.str();
  }
//--------------------------------------------------------------------
  std::pair<string,string> AlexConf::ParseStringPair(const XMLElement* mom, 
                                                    const std::pair<string,string>& tags) 
//--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();

    klog << log4cpp::Priority::DEBUG << " ParseStringPair:: tags.first  = " << tags.first;
    const XMLElement* elem = mom->FirstChildElement (tags.first.c_str()) ;
    string first = elem->GetText();
    klog << log4cpp::Priority::DEBUG << " first = " << first;


    klog << log4cpp::Priority::DEBUG << " ParseStringPair:: tags.second  = " << tags.second;
    const XMLElement*  nextElem = mom->FirstChildElement (tags.second.c_str());
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

    const XMLElement*  nextElem = mom->FirstChildElement (tags.second.c_str());
    //const XMLElement*  nextElem = elem->NextSiblingElement (tags.second.c_str()) ;
    string second = nextElem->GetText();
    klog << log4cpp::Priority::DEBUG << " second = " << second;

    std::pair<int,int> intPair;
    intPair.first = atoi(first.c_str());
    intPair.second = atoi(second.c_str());
    return intPair;

  }

//--------------------------------------------------------------------
  void AlexConf::ParseParamElement(const XMLElement* param) const
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
  void AlexConf::ParseArrayElement(const XMLElement* array) const
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

