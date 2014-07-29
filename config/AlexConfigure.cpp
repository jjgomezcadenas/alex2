
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
using std::pair; 
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
    string rootName = rootElement->Name();

    if (rootName !="AlexConfig")
    {
      klog << log4cpp::Priority::ERROR 
      << "ParseConfiguration::expected root name AlexConfig, found = " 
      << rootName;
      exit (EXIT_FAILURE);
    }

    const XMLElement* algoElement = rootElement->FirstChildElement ("Algos") ;
    fAlgosPathName =  ParseStringPair(algoElement,fStags);

    const XMLElement* dstElement = rootElement->FirstChildElement ("DST") ;
    fDstPathName =  ParseStringPair(dstElement,fStags);

    const XMLElement* histoElement = rootElement->FirstChildElement ("HistoFile") ;
    fHistoPathName =  ParseStringPair(histoElement,fStags);

    const XMLElement* eventElement = rootElement->FirstChildElement ("Events") ;

    std::pair<std::string,std::string> tags;
    tags.first = "runMax";
    tags.second="runDebug";
    fEvents =  ParseIntPair(eventElement,tags);

    const XMLElement* debugElement = rootElement->FirstChildElement ("Debug") ;
  
    const XMLElement* elem = debugElement->FirstChildElement ("level") ;
    fDebug  = elem->GetText();
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

    for (auto i=0; i < fAlgoNames.size(); ++i)
    {
      auto algoPath = fAlgoPath[i];
      auto algoName = fAlgoNames[i];
      vector<DParam> paramVector;

      klog << log4cpp::Priority::DEBUG << " algoPath = " << algoPath;
      fDoc.LoadFile( algoPath.c_str() );
      if (fDoc.ErrorID()!=0) 
      {
        klog << log4cpp::Priority::ERROR 
        << "ParseConfiguration::Failed loading config file error = " 
        << fDoc.ErrorID();
        exit (EXIT_FAILURE);
      }
      XMLElement* rootElement = fDoc.RootElement();  
      const XMLElement* param = rootElement->FirstChildElement ("Param") ;
      
      if (param != NULL)
      {
        vector<DParam> paramVector;

        DParam par = ParseParamElement(param);
        paramVector.push_back(par);
        param = param->NextSiblingElement ("Param") ;
        while (param != NULL)
        {
          DParam par = ParseParamElement(param);
          paramVector.push_back(par);
          param = param->NextSiblingElement ("Param") ;
        }
        fAlgoParam[algoName]=paramVector;
      }

      const XMLElement* array = rootElement->FirstChildElement ("Array") ;

      if (array != NULL)
      {
        vector<DArray> arrayVector;
        DArray ar = ParseArrayElement(array);
        arrayVector.push_back(ar);
        array = array->NextSiblingElement ("Array") ;
        while (array != NULL)
        {
          DArray ar = ParseArrayElement(array);
          arrayVector.push_back(ar);
          array = array->NextSiblingElement ("Array") ;
        }
        fAlgoArray[algoName]=arrayVector;
      }
    
      const XMLElement* h1d = rootElement->FirstChildElement ("H1D") ;

      if (h1d != NULL)
      {
        vector<DH1> h1Vector;
        DH1 h1 = ParseH1DElement(h1d);
        h1Vector.push_back(h1);
        h1d = h1d->NextSiblingElement ("H1D") ;
        while (h1d != NULL)
        {
          DH1 h1 = ParseH1DElement(h1d);
          h1Vector.push_back(h1);
          h1d = h1d->NextSiblingElement ("H1D") ;
        }
        fAlgoH1D[algoName]=h1Vector;
      }

      const XMLElement* h2d = rootElement->FirstChildElement ("H2D") ;

      if (h2d != NULL)
      {
        vector<DH2> h2Vector;
        DH2 h2 = ParseH2DElement(h2d);
        h2Vector.push_back(h2);
        h2d = h2d->NextSiblingElement ("H2D") ;
        while (h2d != NULL)
        {
          vector<DH2> h2Vector;
          DH2 h2 = ParseH2DElement(h2d);
          h2Vector.push_back(h2);
          h2d = h2d->NextSiblingElement ("H2D") ;
        }
        fAlgoH2D[algoName]=h2Vector;
      }
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
    string rootName = rootElement->Name();

    if (rootName !="AlgoConfig")
    {
      klog << log4cpp::Priority::ERROR 
      << "ParseConfiguration::expected root name AlexConfig, found = " 
      << rootName;
      exit (EXIT_FAILURE);
    }

    const XMLElement* algoElement = rootElement->FirstChildElement ("Algo") ;
    
    if (algoElement != NULL)
    {
      pair<string,string> algoPathName = ParseStringPair(algoElement,fStags);
      string algoPath = PathFromStrings(algoPathName.first,algoPathName.second);
      klog << log4cpp::Priority::DEBUG << " algoPath " << algoPath;
      fAlgoNames.push_back(algoPathName.second);
      algoPath = MergeStrings(algoPath,".xml");
      fAlgoPath.push_back(algoPath);

      const XMLElement*  nextAlgo = algoElement->NextSiblingElement ("Algo") ;

      while (nextAlgo != NULL)
      {
        algoPathName = ParseStringPair(nextAlgo,fStags);
        algoPath = PathFromStrings(algoPathName.first,algoPathName.second);
        algoPath = MergeStrings(algoPath,".xml");
        fAlgoNames.push_back(algoPathName.second);
        fAlgoPath.push_back(algoPath);
        nextAlgo = nextAlgo->NextSiblingElement ("Algo") ;
      }
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
  std::string AlexConf::SerializeAlgoParam() const
//--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    
    ostringstream s;
    s << std::endl;
    
    for (auto& kv : fAlgoParam) 
    {
      string algoName = kv.first;
      klog << log4cpp::Priority::DEBUG << " algo name " << algoName;
      s << "Algo: =" << algoName << endl;
      
      std::vector<alex::DParam> vparam = kv.second;
      for (auto param : vparam)
      {
        klog << log4cpp::Priority::DEBUG << " ++Param++ ";
        klog << log4cpp::Priority::DEBUG <<param.Serialize();

        s << param.Serialize() << endl; 
      }
    }
    //s<< std::ends;
    return s.str();
  }

//--------------------------------------------------------------------
  std::string AlexConf::SerializeAlgoArray() const
//--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    
    ostringstream s;
    s << std::endl;
    
    for (auto& kv : fAlgoArray) 
    {
      string algoName = kv.first;
      klog << log4cpp::Priority::DEBUG << " algo name " << algoName;
      s << "Algo: =" << algoName << endl;
      
      std::vector<alex::DArray> vparam = kv.second;
      for (auto param : vparam)
      {
        klog << log4cpp::Priority::DEBUG << " ++Array++ ";
        klog << log4cpp::Priority::DEBUG <<param.Serialize();

        s << param.Serialize() << endl; 
      }
    }
    //s<< std::ends;
    return s.str();
  }

//--------------------------------------------------------------------
  std::string AlexConf::SerializeAlgoH1D() const
//--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    
    ostringstream s;
    s << std::endl;
    
    for (auto& kv : fAlgoH1D) 
    {
      string algoName = kv.first;
      klog << log4cpp::Priority::DEBUG << " algo name " << algoName;
      s << "Algo: =" << algoName << endl;
      
      std::vector<alex::DH1> vparam = kv.second;
      for (auto param : vparam)
      {
        klog << log4cpp::Priority::DEBUG << " ++H1D++ ";
        klog << log4cpp::Priority::DEBUG <<param.Serialize();

        s << param.Serialize() << endl; 
      }
    }
    //s<< std::ends;
    return s.str();
  }
//--------------------------------------------------------------------
  std::string AlexConf::SerializeAlgoH2D() const
//--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    
    ostringstream s;
    s << std::endl;
    
    for (auto& kv : fAlgoH2D) 
    {
      string algoName = kv.first;
      klog << log4cpp::Priority::DEBUG << " algo name " << algoName;
      s << "Algo: =" << algoName << endl;
      
      std::vector<alex::DH2> vparam = kv.second;
      for (auto param : vparam)
      {
        klog << log4cpp::Priority::DEBUG << " ++H2D++ ";
        klog << log4cpp::Priority::DEBUG <<param.Serialize();

        s << param.Serialize() << endl; 
      }
    }
    //s<< std::ends;
    return s.str();
  }
//--------------------------------------------------------------------
  std::string AlexConf::WriteAConfHeader() const
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
    //s << std::ends;
    return s.str();
  }

//--------------------------------------------------------------------
  std::string AlexConf::WriteAConfCPP() const
//--------------------------------------------------------------------
  {
    ostringstream s;
    s<<"\n// Generated by AlexConf: do not edit" << endl;

    s<<"#include "<< '"'<<"AConf.hh"<< '"'<< endl;
    s<<"using namespace alex ;" << endl;
    s<<"AAConf::AAConf()" << endl;
    s<<"{" << endl;
    s<<"  fAlgoPathName.first=" << '"'<<fAlgosPathName.first << '"' <<";" <<endl;
    s<<"  fAlgoPathName.second=" << '"'<< fAlgosPathName.second << '"' <<";" <<endl;
    s<<"  fDstPathName.first=" << '"'<<fDstPathName.first << '"' <<";" <<endl;
    s<<"  fDstPathName.second=" << '"'<< fDstPathName.second << '"' <<";" <<endl;
    s<<"  fHistoPathName.first=" << '"'<< fHistoPathName.first << '"' <<";" <<endl;
    s<<"  fHistoPathName.second=" << '"'<< fHistoPathName.second << '"' <<";" <<endl;
    s<<"  fEvents.first=" << fEvents.first <<";" <<endl;
    s<<"  fEvents.second=" << fEvents.second <<";" <<endl;
    s<<"  fDebug=" << '"'<< fDebug << '"' <<";" <<endl;
    s<<"}" << endl;
    
    //s << std::ends;
    return s.str();
  }
//--------------------------------------------------------------------
  std::string AlexConf::WriteAlgoHeader()
//--------------------------------------------------------------------
  {
    ostringstream s;

    s<<"\n#ifndef ALGOHEAD_" << endl;
    s<<"#define ALGOHEAD_" << endl;
    s<<"// Generated by AlexConf: do not edit" << endl;

    s<<"#include <string>" << endl;
    s<<"#include <vector>" << endl;
    s<<"#include <utility>" << endl;
    s<<"#include <memory>" << endl;
    s<<"#include <map>" << endl;
    s<<"#include <TF1.h>" << endl;
    s<<"#include <TF2.h>" << endl;
    s<<"#include <alex/IAlgorithm.h>" << endl;

    s<<"namespace alex {" << endl;

    
    for (auto algoName : fAlgoNames) 
    {
      
      s << "  class " << algoName << ": public IAlgorithm {" << endl;
      s << "  public:" << endl;
      s << "    bool Init() ;"<< endl;
      s << "    bool Execute() ;"<< endl;
      s << "    bool End() ;"<< endl;
      s << "    std::string  Name() const {return fName;}"<< endl;
      s << "    void SetName(std::string name) {fName = name;}"<< endl;
      s << "  private:"<< endl;
      std::vector<alex::DParam> vparam = fAlgoParam[algoName];
      for (auto param : vparam)
      {
        if(param.DataType()=="double")
          s << "    double "<< param.Name() << ";" <<endl;
        else if(param.DataType()=="int")
          s << "    int "<< param.Name() << ";" <<endl;
        else
          s << "    std::string "<< param.Name() << ";" <<endl;
      }

      std::vector<alex::DArray> varray = fAlgoArray[algoName];
      for (auto array : varray)
      {
        if(array.DataType()=="double")
          s << "    std::vector<double> "<< array.Name() << ";" <<endl;
        else if(array.DataType()=="int")
          s << "    std::vector<int> "<< array.Name() << ";" <<endl;
        else
          s << "    std::vector<std::string> "<< array.Name() << ";" <<endl;
      }
      std::vector<alex::DH1> vh1 = fAlgoH1D[algoName];
      for (auto h1 : vh1)
      {
        s << "    TF1* "<< h1.Name() << ";" <<endl;
      }

      std::vector<alex::DH2> vh2 = fAlgoH2D[algoName];
      for (auto h2 : vh2)
      {
        s << "    TF2* "<< h2.Name() << ";" <<endl;
      }

      s << "  };"<< endl;
    }
    s << "}"<< endl;
    s << "#endif "<< endl;
    
    //s << std::ends;
    return s.str();
  }

//--------------------------------------------------------------------
  vector<string> AlexConf::WriteAlgoHeaders()
//--------------------------------------------------------------------
  {
    vector<string> algoHeaders;
    for (auto algoName : fAlgoNames) 
    {
      ostringstream s;

      s<<"\n#ifndef ALGO"<<algoName<<"_" << endl;
      s<<"#define ALGO"<<algoName<<"_" << endl;
      s<<"// Generated by AlexConf: do not edit" << endl;

      s<<"#include <string>" << endl;
      s<<"#include <vector>" << endl;
      s<<"#include <utility>" << endl;
      s<<"#include <memory>" << endl;
      s<<"#include <map>" << endl;
      s<<"#include <TH1F.h>" << endl;
      s<<"#include <TH2F.h>" << endl;
      s<<"#include <alex/IAlgorithm.h>" << endl;

      s<<"namespace alex {" << endl;
      s << "  class " << algoName << ": public IAlgorithm {" << endl;
      s << "  public:" << endl;
      s << "    "<< algoName << "();" << endl;
      s << "    ~"<< algoName << "(){}" << endl;
      s << "    bool Init() ;"<< endl;
      s << "    bool Execute() ;"<< endl;
      s << "    bool End() ;"<< endl;
      s << "    std::string  Name() const {return fName;}"<< endl;
      s << "    void SetName(std::string name) {fName = name;}"<< endl;
      s << "  private:"<< endl;
      std::vector<alex::DParam> vparam = fAlgoParam[algoName];
      for (auto param : vparam)
      {
        if(param.DataType()=="double")
          s << "    double "<< param.Name() << ";" <<endl;
        else if(param.DataType()=="int")
          s << "    int "<< param.Name() << ";" <<endl;
        else
          s << "    std::string "<< param.Name() <<";" <<endl;
      }

      std::vector<alex::DArray> varray = fAlgoArray[algoName];
      for (auto array : varray)
      {
        if(array.DataType()=="double")
          s << "    std::vector<double> "<< array.Name() << ";" <<endl;
        else if(array.DataType()=="int")
          s << "    std::vector<int> "<< array.Name() << ";" <<endl;
        else
          s << "    std::vector<std::string> "<< array.Name() << ";" <<endl;
      }
      std::vector<alex::DH1> vh1 = fAlgoH1D[algoName];
      for (auto h1 : vh1)
      {
        s << "    TH1F* "<< h1.Name() << ";" <<endl;
      }

      std::vector<alex::DH2> vh2 = fAlgoH2D[algoName];
      for (auto h2 : vh2)
      {
        s << "    TH2F* "<< h2.Name() << ";" <<endl;
      }

      s << "  };"<< endl;
      s << "}"<< endl;
      s << "#endif "<< endl;
      algoHeaders.push_back(s.str());
    }
    return algoHeaders;
  }

//--------------------------------------------------------------------
  std::string AlexConf::WriteAlgoCPP()
//--------------------------------------------------------------------
  {
    ostringstream s;

    s<<"// Generated by AlexConf: do not edit" << endl;
    for (auto algoName : fAlgoNames) 
    {
      s<<"#include "<< '"'<<algoName<<".hh"<< '"'<< endl;
    }

    s<<"namespace alex {" << endl;

    for (auto algoName : fAlgoNames) 
    {
      s << "\n\n  "<<algoName<<"::"<<algoName<<"()"<<endl;
      s << "  {" << endl;
      
      std::vector<alex::DParam> vparam = fAlgoParam[algoName];
      for (auto param : vparam)
        if(param.DataType()=="string")
          s << "   "<< param.Name() << " = " <<'"' <<param.Value() <<'"'<< ";" << endl;
        else
          s << "   "<< param.Name() << " = " <<param.Value() << ";" << endl;
          
      std::vector<alex::DArray> varray = fAlgoArray[algoName];
      for (auto array : varray)
      {
        std::vector<std::string> values = array.Split() ;
        s << "   "<< array.Name() << ".reserve(" <<array.Dim() << ");" << endl;
        for (auto value: values)
          s << "   "<< array.Name() << ".push_back(" <<value << ");" << endl;
      }

      std::vector<alex::DH1> vh1 = fAlgoH1D[algoName];
      for (auto h1 : vh1)
      {
        s << "    "<< h1.Name() << "= new TH1F("<<'"'<<h1.Name() <<'"' <<"," <<
          '"'<<h1.Title() <<'"' <<"," <<
          h1.Nbinsx() <<"," <<
          h1.Xlow() <<"," <<
          h1.Xup() << ");" << endl;
      }

      std::vector<alex::DH2> vh2 = fAlgoH2D[algoName];
      for (auto h2 : vh2)
      {
        s << "    "<< h2.Name() << "= new TH2F("<<'"'<<h2.Name() <<'"' <<"," <<
          '"'<<h2.Title() <<'"' <<"," <<
          h2.Nbinsx() <<"," <<
          h2.Xlow() <<"," <<
          h2.Xup() <<"," <<
          h2.Nbinsy() <<"," <<
          h2.Ylow() <<"," <<
          h2.Yup()  << ");" << endl;
      }
      
      s << "  }" << endl;
    }
    s << "}" << endl;
    //s << std::ends;
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
  DParam AlexConf::ParseParamElement(const XMLElement* param) const
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

    auto par = DParam();
    par.SetData(textNameParam,dataTypeParamText,ValueParamText);
    return par;
  }
//--------------------------------------------------------------------
  DArray  AlexConf::ParseArrayElement(const XMLElement* array) const
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

    auto ar = DArray();
    ar.SetData(textName,textDataType,textDim,textValue);
    return ar;
  }

//--------------------------------------------------------------------
  DH1 AlexConf::ParseH1DElement(const XMLElement* h1d) const
//--------------------------------------------------------------------
  {
    //   <H1D>
    //     <name>fHY</name>
    //     <title>Y distribution</title>
    //     <nbinsx>50</nbinsx>
    //     <xlow>-10.0</xlow>
    //     <xup>10.0</xup>
    //   </H1D>
    
    log4cpp::Category& klog = log4cpp::Category::getRoot();

    const XMLElement* name = h1d->FirstChildElement ("name") ;
    string textName = name->GetText();
    klog << log4cpp::Priority::DEBUG << " H1D name text = " << textName;

    const XMLElement* title = h1d->FirstChildElement ("title") ;
    string textTitle= title->GetText();
    klog << log4cpp::Priority::DEBUG << "H1D title text =" << textTitle;

    const XMLElement* nbinsx = h1d->FirstChildElement ("nbinsx") ;
    string textNbinsx= nbinsx->GetText();
    klog << log4cpp::Priority::DEBUG << "nbinsx text =" << textNbinsx;

    const XMLElement* xlow = h1d->FirstChildElement ("xlow") ;
    string textXlow= xlow->GetText();
    klog << log4cpp::Priority::DEBUG << "xlow text =" << textXlow;

    const XMLElement* xup = h1d->FirstChildElement ("xup") ;
    string textXup= xup->GetText();
    klog << log4cpp::Priority::DEBUG << "xup text =" << textXup;

    auto h1 = DH1();
    h1.SetData(textName,textTitle,textNbinsx,textXlow,textXup);
    return h1;
  }
//--------------------------------------------------------------------
  DH2 AlexConf::ParseH2DElement(const XMLElement* h2d) const
//--------------------------------------------------------------------
  {
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
    log4cpp::Category& klog = log4cpp::Category::getRoot();

    const XMLElement* name = h2d->FirstChildElement ("name") ;
    string textName = name->GetText();
    klog << log4cpp::Priority::DEBUG << " H2D name text = " << textName;

    const XMLElement* title = h2d->FirstChildElement ("title") ;
    string textTitle= title->GetText();
    klog << log4cpp::Priority::DEBUG << "H2D title text =" << textTitle;

    const XMLElement* nbinsx = h2d->FirstChildElement ("nbinsx") ;
    string textNbinsx= nbinsx->GetText();
    klog << log4cpp::Priority::DEBUG << "nbinsx text =" << textNbinsx;

    const XMLElement* nbinsy = h2d->FirstChildElement ("nbinsy") ;
    string textNbinsy= nbinsy->GetText();
    klog << log4cpp::Priority::DEBUG << "nbinsy text =" << textNbinsy;

    const XMLElement* xlow = h2d->FirstChildElement ("xlow") ;
    string textXlow= xlow->GetText();
    klog << log4cpp::Priority::DEBUG << "xlow text =" << textXlow;

    const XMLElement* xup = h2d->FirstChildElement ("xup") ;
    string textXup= xup->GetText();
    klog << log4cpp::Priority::DEBUG << "xup text =" << textXup;

    const XMLElement* ylow = h2d->FirstChildElement ("ylow") ;
    string textYlow= ylow->GetText();
    klog << log4cpp::Priority::DEBUG << "ylow text =" << textYlow;

    const XMLElement* yup = h2d->FirstChildElement ("yup") ;
    string textYup= yup->GetText();
    klog << log4cpp::Priority::DEBUG << "yup text =" << textYup;

    auto h2 = DH2();
    h2.SetData(textName,textTitle,textNbinsx,textXlow,textXup,
      textNbinsy,textYlow,textYup);

    return h2;

  }
}

