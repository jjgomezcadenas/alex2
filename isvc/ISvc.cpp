
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
    for (auto it = 0; it < (int) fIevt->Tracks().size(); ++it)
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

//
// maxDist: the maximum distance a hit can be from some hit in a track
//  to still be considered part of that track
//--------------------------------------------------------------------
  void IreneManager::CreateTracks(double maxDist)
//--------------------------------------------------------------------
  {
    // Clear the current list of IBeta objects.
    fIBetas.clear();

    // Clear the list of minimum distances.
    minDistList.clear();

    // Get the hits corresponding to the tracks produced by the electrons with the most
    //  and second-most momentum. 
    std::pair<IHits, IHits> electronhits = ISvc::Instance().GetPMaxElectronsHits();
 
    // Get all hits.
    IHits truehits = ISvc::Instance().GetTrueHits();
    
    // Create the IBeta object.
    IBeta * current_track = new IBeta();

    // Add all the hits from the most and second-most energetic
    //  electrons as the seed.
    for(int i = 0; i < electronhits.first.size(); i++) current_track->AddHit(electronhits.first[i]);
    for(int i = 0; i < electronhits.second.size(); i++) current_track->AddHit(electronhits.second[i]);

    // First loop over all hits and calculate the minimum distance to the primary track.
    for(int i = 0; i < truehits.size(); i++)
    {
      IHit ih = truehits[i];

      // Compute the minimum distance between all the hits in the
      // track thus far and the current hit.
      double mindist = ComputeMinDist(current_track->GetHit(), ih);
      
      // Record the minimum distance if the hit is not in the primary track.
      if(mindist > 1.0e-8)
      {
        minDistList.push_back(mindist);
      }
      
    }

    // Loop over the hits, keeping track of the number of hits that
    //  were not added to the current track (variable i).
    int i = 0;
    while(truehits.size() != 0)
    {
      //std::cout << "-- Looping with truehits.size = " << truehits.size() 
      //          << " and i = " << i << std::endl;

      IHit ih = truehits[i];

      // Compute the minimum distance between all the hits in the
      // track thus far and the current hit.
      double mindist = ComputeMinDist(current_track->GetHit(), ih);
      //std::cout << "Found mindist = " << mindist << std::endl;
      
      // If mindist is zero within some tolerance, we have probably found the same hit: remove it.
      if(mindist < 1.0e-8)
      {
        truehits.erase(truehits.begin()+i);
      }
      // Otherwise, if mindist is less than the chosen track separation, add to the current track.
      else if(mindist < maxDist)
      {
  
        // This hit lies in the curent track; add it.
        current_track->AddHit(ih);
        //std::cout << "Added close hit at (X, Y, Z) = (" << ih.first.X() << ", " << ih.first.Y() 
        //      << ", " << ih.first.Z() << ")" << std::endl;

        // Remove from the list of all hits.
        truehits.erase(truehits.begin()+i);

        // Reset the hit counter.
        i = 0;
      }
      else
      {
        // This hit lies outside the current track; advance to the next one.
        i++;
      }

      // If all remaining hits lie outside the current track, record
      //  the current track object, and create a new current track.
      if(i >= truehits.size() && truehits.size() != 0)
      {

        // Save the current track.
        fIBetas.push_back(current_track);

        // Start the next track with the first hit in the list of remaining hits.
        current_track = new IBeta();
        current_track->AddHit(truehits[0]);
        truehits.erase(truehits.begin());
        i = 0;
      }
    }

    // Push back the final track.
    fIBetas.push_back(current_track);
  }

/*
// 
// Create IBeta objects from hits without the use of the electron track as seed.
//--------------------------------------------------------------------
  void IreneManager::CreateUnseededTracks()
//--------------------------------------------------------------------
  {

    // Clear the current list of IBeta objects.
    fIBetas.clear();

    // Get all hits.
    IHits truehits =ISvc::Instance().GetTrueHits();

    // Create the IBeta object.
    IBeta * current_track = new IBeta();

    // Add the first hit as the seed.
    current_track->AddHit(truehits[0]);
    truehits.erase(truehits.begin());
    //std::cout << "Added hit at (X, Y, Z) = (" << truehits[0].first.X() << ", " << truehits[0].first.Y() 
    //          << ", " << truehits[0].first.Z() << ")" << std::endl;

    // Loop over the hits, keeping track of the number of hits that
    //  were not added to the current track (variable i).
    int i = 0;
    while(truehits.size() != 0)
    {
      //std::cout << "-- Looping with truehits.size = " << truehits.size() 
      //          << " and i = " << i << std::endl;

      IHit ih = truehits[i];

      // Compute the minimum distance between all the hits in the
      // track thus far and the current hit.
      double mindist = ComputeMinDist(current_track->GetHit(), ih);
      //std::cout << "Found mindist = " << mindist << std::endl;
      if(mindist < fMaxDist)
      {
        // This hit lies in the curent track; add it.
        current_track->AddHit(ih);
        //std::cout << "Added close hit at (X, Y, Z) = (" << ih.first.X() << ", " << ih.first.Y() 
        //      << ", " << ih.first.Z() << ")" << std::endl;

        // Remove from the list of all hits.
        truehits.erase(truehits.begin());

        // Reset the hit counter.
        i = 0;
      }
      else
      {
        // This hit lies outside the current track; advance to the next one.
        i++;
      }

      // If all remaining hits lie outside the current track, record and
      //  close the current track object.
      if(i >= truehits.size() && truehits.size() != 0)
      {
        // Save the current track.
        fIBetas.push_back(current_track);

        // Start the next track with the first hit in the list of remaining hits.
        current_track->Clear();
        current_track->AddHit(truehits[0]);
        truehits.erase(truehits.begin());
        i = 0;
      }
    }
  }
*/

//--------------------------------------------------------------------
  double IreneManager::ComputeMinDist(IHits ihs1, IHits ihs2)
//--------------------------------------------------------------------
  {
    double mindist = 1.e9;

    for (int i = 0; i<ihs1.size(); i++)
    {
      for (int j = 0; j<ihs2.size(); j++)
      {
        double diffX = ihs1[i].first.X()-ihs2[j].first.X();
        double diffY = ihs1[i].first.Y()-ihs2[j].first.Y();
        double diffZ = ihs1[i].first.Z()-ihs2[j].first.Z();
        double dist = sqrt(diffX*diffX+diffY*diffY+diffZ*diffZ);
       
        if(dist<mindist)
          mindist = dist;
      }
    }

    return mindist;
  }

//--------------------------------------------------------------------
  double IreneManager::ComputeMinDist(std::vector<const alex::Hit*> ihs1, IHit ih)
//--------------------------------------------------------------------
  {
    double mindist = 1.e9;

    for (int i = 0; i<ihs1.size(); i++)
    {
      double diffX = ihs1[i]->XYZ().X()-ih.first.X();
      double diffY = ihs1[i]->XYZ().Y()-ih.first.Y();
      double diffZ = ihs1[i]->XYZ().Z()-ih.first.Z();
      double dist = sqrt(diffX*diffX+diffY*diffY+diffZ*diffZ);
      //std::cout << "Comparing hit (" << ih.first.X() << ", " << ih.first.Y()
      //        << ", " << ih.first.Z() << ") with hit at (X, Y, Z) = (" << ihs1[i]->XYZ().X() 
      //        << ", " << ihs1[i]->XYZ().Y() << ", " << ihs1[i]->XYZ().Z() << "): dist = " << dist
      //        << std::endl;


      if(dist<mindist)
        mindist = dist;
    }

    return mindist;
  }

}
