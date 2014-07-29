
// ----------------------------------------------------------------------------
//  $Id: Alex.cpp 
//
//  Author:  <gomez@mail.cern.ch>
//  Created: July 2014
//  
//  Copyright (c) 2014 NEXT Collaboration
// ---------------------------------------------------------------------------- 

#include <alex/ISvc.h>
#include <alex/LogUtil.h>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <cstdio>
#include <cstring>

using std::string; 
using std::cout; 
using std::endl; 
using std::ostream;
using std::vector;

namespace alex {

//--------------------------------------------------------------------
  void IreneManager::InitDstFile(std::string fileName)
//--------------------------------------------------------------------
  {
    fIfile = new TFile(fileName.c_str(), "READ");
    fEvtTree = dynamic_cast<TTree*>(fIfile->Get("EVENT"));
  }
//--------------------------------------------------------------------
  int IreneManager::DstEntries()
//--------------------------------------------------------------------
  {
    return (int) fEvtTree->GetEntries();
  }

//--------------------------------------------------------------------
  int IreneManager::DstGetEntry(int ivt) 
//--------------------------------------------------------------------
  {
    return fEvtTree->GetEntry(ivt);
  }
//--------------------------------------------------------------------
  void IreneManager::Init(string debugLevel)
//--------------------------------------------------------------------
  {
    SetDebugLevel(debugLevel);

  }
//--------------------------------------------------------------------
  void IreneManager::LoadEvent(const irene::Event* ievt)
//--------------------------------------------------------------------
  {
    fIevt = ievt;
  }
//--------------------------------------------------------------------
  const irene::Event& IreneManager::GetEvent()
//--------------------------------------------------------------------
  {
    return *fIevt;
  }
}

