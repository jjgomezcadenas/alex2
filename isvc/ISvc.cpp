
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
  void IreneManager::Init(string debugLevel)
//--------------------------------------------------------------------
  {
    SetDebugLevel(debugLevel);
  }
//--------------------------------------------------------------------
  void IreneManager::InitDst(std::string fileName, const irene::Event* ievt)
//--------------------------------------------------------------------
  {
    fIfile = new TFile(fileName.c_str(), "READ");
    fEvtTree = dynamic_cast<TTree*>(fIfile->Get("EVENT"));
    fEvtTree->SetBranchAddress("EventBranch", &ievt);
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
  void IreneManager::LoadEvent(const irene::Event* ievt)
//--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::DEBUG << " IreneManager::LoadEvent =" ;

    fIevt = ievt;

    klog << log4cpp::Priority::DEBUG << " IreneManager:: call FetchElectrons()" ;
    FetchElectrons();
    FetchPMaxElectrons();
    fTrueHits.clear();
    fIevt->FillHitVector(fTrueHits, "ACTIVE");
    fIreneTracks = fIevt->Tracks();

    for (auto track: fIreneTracks)
    {
      auto particle = track->GetParticle();
      if (particle->GetParticleID() == fBetasMax.first->GetParticleID())
      {
        fBetasMaxHits.first = track->GetHits();
      }
      else if (particle->GetParticleID() == fBetasMax.second->GetParticleID())
      {

        fBetasMaxHits.second = track->GetHits();
      }
    }
  }
//--------------------------------------------------------------------
  const irene::Event& IreneManager::GetEvent()
//--------------------------------------------------------------------
  {
    return *fIevt;
  }

//--------------------------------------------------------------------
  void IreneManager::FetchElectrons() 
//--------------------------------------------------------------------
  {
    fElectrons.clear();
    fBetas.clear();
    for (auto it = 0; it < fIevt->Tracks().size(); ++it)
    {  
      const irene::Particle* ip= fIevt->Tracks().at(it)->GetParticle();
    
      //if (ip->GetPDGcode() == 11 and ip->IsPrimary() == true) 
      if (ip->GetPDGcode() == 11) 
      {
        fElectrons.push_back(ip);
        if (ip->IsPrimary() == true)
          fBetas.push_back(ip);
      }
    } 
  }
//--------------------------------------------------------------------
  void IreneManager::FetchPMaxElectrons() 
//--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::DEBUG << " IreneManager::GetPMaxElectrons()" ;

    
    double pmax = 0;
    int imax=-1;
    int i=0;
    for(auto beta : fBetas)
    {
      if (beta->Momentum() > pmax)
      {
        pmax = beta->Momentum();
        imax=i;
      }
      i++;
    }

    klog << log4cpp::Priority::DEBUG << " imax =" << imax << " pmaxx = " << pmax;

    fBetasMax.first = fBetas.at(imax);
    if (GetNumberOfPrimaryElectrons() <2 )
      fBetasMax.second = NULL;
    else
    {
      pmax = 0;
      int imax2=-1;
      i=-1;
      for(auto beta : fBetas)
      {
        i++;
        if (i==imax) continue;
        if (beta->Momentum() > pmax)
        {
          pmax = beta->Momentum();
          imax2=i;
        }
      }
      klog << log4cpp::Priority::DEBUG << " imax2 =" << imax2 << " pmax = " << pmax;
      fBetasMax.second =fBetas.at(imax2);
    }
  }
}

