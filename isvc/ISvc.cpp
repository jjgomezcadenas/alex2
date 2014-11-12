
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
    fRandom = new TRandom2();
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
      else if (GetNumberOfPrimaryElectrons() >= 2 && particle->GetParticleID() == fBetasMax.second->GetParticleID())
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

    // Checksum for all hits.
    double echeck = 0.;
    for(int i = 0; i < (int) truehits.size(); i++) {
      echeck += truehits[i].second;
    } 
    std::cout << "TOTAL energy in all hits is = " << echeck << std::endl;   
 
    // Create the IBeta object for the primary track.
    IBeta * current_track = new IBeta();

    // Add all the hits from the most and second-most energetic
    //  electrons as the seed.
    echeck = 0.;
    for(int i = (int) (electronhits.first.size()-1); i > 0; i--) {
      current_track->AddHit(electronhits.first[i]);
      echeck += electronhits.first[i].second;
    }
    if(GetNumberOfPrimaryElectrons() >= 2) {
      for(int i = 0; i < (int) electronhits.second.size(); i++) {
        current_track->AddHit(electronhits.second[i]);
        echeck += electronhits.second[i].second;
      }
    }
    std::cout << "SUM of deltaE in most and second-most energetic electrons = " << echeck << std::endl;

    // Find the blob centers.
    IHits trkFirst = electronhits.first;
    IHits trkSecond = electronhits.second;
    TLorentzVector cblob1, cblob2;
    cblob1 = trkFirst[trkFirst.size()-1].first;
    if(ISvc::Instance().GetNumberOfPrimaryElectrons() >= 2) cblob2 = trkSecond[trkSecond.size()-1].first;
    else cblob2 = trkFirst[0].first;

    // Set the vertices accordingly.
    TVector3 vt1(0.,0.,0.);
    TVector3 vt2(0.,0.,0.);
    // If we have a double-beta event, set the V0 to the beginning of the highest-momentum track,
    //  and the ends to the first and second blobs.
    if(ISvc::Instance().GetNumberOfPrimaryElectrons() >= 2) {
      vt1.SetX(trkFirst[0].first.X());
      vt1.SetY(trkFirst[0].first.Y());
      vt1.SetZ(trkFirst[0].first.Z());
      current_track->SetV0(vt1);
    }
    vt1.SetX(cblob1.X()); vt2.SetX(cblob2.X());
    vt1.SetY(cblob1.Y()); vt2.SetY(cblob2.Y());
    vt1.SetZ(cblob1.Z()); vt2.SetZ(cblob2.Z());
    current_track->SetVX(vt1,vt2);

    // Loop over all hits and calculate the minimum distance to the primary track.
    for(int i = 0; i < (int) truehits.size(); i++)
    {
      IHit ih = truehits[i];

      // Compute the minimum distance between all the hits in the
      // track thus far and the current hit.
      int imin, iside;
      double mindist = ComputeMinDist(current_track->GetHit(), ih, imin, iside);
      
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
      int imin, iside;
      double mindist = ComputeMinDist(current_track->GetHit(), ih, imin, iside);
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
        Hit trackHit(ih);
        trackHit.SetMinDist(mindist);  // save the minimum distance so, for the primary track, it does not
                                       // need to be recomputed later
        current_track->InsertHit(trackHit,std::max(imin,iside));
        //std::cout << "Added hit to IBeta " << fIBetas.size() << " with mindist = " << mindist << std::endl;
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
      if(i >= (int) truehits.size() && truehits.size() != 0)
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

//
// coreDist: the maximum distance (in mm) a hit can be from the main track and
//  still be considered part of the core
// zPlaneWidth: the z-plane width (in mm)
//--------------------------------------------------------------------
  void IreneManager::CreateRTracks(double coreDist, int sparseWidth, double blobRadius)
//--------------------------------------------------------------------
  {
    double echeck = 0.;

    // Get the hits corresponding to the tracks produced by the electrons with the most
    //  and second-most momentum. 
    std::pair<IHits, IHits> electronhits = ISvc::Instance().GetPMaxElectronsHits();
 
    // Get all hits in the main IBeta object.
    std::vector<const Hit*> trackHits = fIBetas[0]->GetHit();
    IHits photonHits;  // to temporarily keep track of photon hits

    // -----------------------------------------------------------------------------
    // Create the RBeta object
    // -----------------------------------------------------------------------------

    // Create the RBeta object and copy all hits in the main track.
    fRBeta = new RBeta();
    //double zminCore = 0., zmaxCore = 0.;    // save the minimum and maximum core z-coordinate for later use
    echeck = 0.;
    for(int i = 0; i < (int) trackHits.size(); i++) {

      // Add this hit to the RBeta.
      fRBeta->AddHit(*(trackHits[i]));

      // If it is a core hit, add it to the list of core hits.
      if(trackHits[i]->MinDist() < coreDist) {
        //double zcore = trackHits[i]->XYZ().Z();
        //if(i == 0 || zcore < zminCore) zminCore = zcore;
        //if(i == 0 || zcore > zmaxCore) zmaxCore = zcore;
        fRBeta->AddCoreHit(*(trackHits[i]));
      }
      // If not, add it to a temporary list of photon hits.
      else {
        IHit photonHit;
        photonHit.first.SetXYZT(trackHits[i]->XYZ().X(),trackHits[i]->XYZ().Y(),trackHits[i]->XYZ().Z(),0.);
        photonHit.second = trackHits[i]->Edep();
        photonHits.push_back(photonHit);
      }
      echeck += trackHits[i]->Edep();
    }
    std::cout << "Total trackHits energy = " << echeck << std::endl;

    // -------------------------------------------------------------------------------------------------------
    // Create the photon tracks from the photon hits.
    // -------------------------------------------------------------------------------------------------------
    
    if(photonHits.size() > 0) {

      // Create an IBeta object to store the current photon track.
      IBeta * current_track = new IBeta();
      current_track->AddHit(photonHits[0]);
      photonHits.erase(photonHits.begin());
      double ephot = 0.;
      double xavg = 0., yavg = 0., zavg = 0.;

      // Loop over the hits, keeping track of the number of hits that
      //  were not added to the current track (variable i).
      int i = 0;
      while(photonHits.size() != 0)
      {
        //std::cout << "-- Looping with truehits.size = " << truehits.size() 
        //          << " and i = " << i << std::endl;
        IHit ih = photonHits[i];

        // Compute the minimum distance between all the hits in the
        // track thus far and the current hit.
        int imin, iside;
        double mindist = ComputeMinDist(current_track->GetHit(), ih, imin, iside);
        //std::cout << "Found mindist = " << mindist << std::endl;
      
        // If mindist is less than the chosen track separation, add to the current track.
        if(mindist < coreDist)
        {
  
          // This hit lies in the curent track; add it.
          Hit trackHit(ih);
          trackHit.SetMinDist(mindist);  // save the minimum distance so, for the primary track, it does not
                                         // need to be recomputed later
          current_track->AddHit(trackHit);
          xavg += trackHit.XYZ().X()*trackHit.Edep();
          yavg += trackHit.XYZ().Y()*trackHit.Edep();
          zavg += trackHit.XYZ().Z()*trackHit.Edep();
          ephot += trackHit.Edep();
          //std::cout << "Added close hit at (X, Y, Z) = (" << ih.first.X() << ", " << ih.first.Y() 
          //      << ", " << ih.first.Z() << ")" << std::endl;

          // Remove from the list of all hits.
          photonHits.erase(photonHits.begin()+i);

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
        if(i >= (int) photonHits.size() && photonHits.size() != 0)
        {

          // Save the current track.
          xavg /= ephot; yavg /= ephot; zavg /= ephot;
          current_track->CreateBlob(xavg,yavg,zavg,ephot);
          fRBeta->AddPhoton(*current_track);
          xavg = 0; yavg = 0; zavg = 0; ephot = 0;

          // Start the next track with the first hit in the list of remaining hits.
          current_track = new IBeta();
          current_track->AddHit(photonHits[0]);
          photonHits.erase(photonHits.begin());
          i = 0;
        }
      }

      // Push back the final photon track.
      xavg /= ephot; yavg /= ephot; zavg /= ephot;
      current_track->CreateBlob(xavg,yavg,zavg,ephot);
      fRBeta->AddPhoton(*current_track);
      xavg = 0; yavg = 0; zavg = 0; ephot = 0;
    }

    // -----------------------------------------------------------------------------
    // Create the effective hits 
    // -----------------------------------------------------------------------------

    // Sparse the hits.
    int ihit = 0; echeck = 0.;
    std::vector<const Hit*> coreHits = fRBeta->GetCoreHits(); 
    while(ihit < (int) coreHits.size()) {

      // Keep track of the sparsing counter.
      int sp = ihit;

      // Create an effective (sparsed) hit.
      double xeff = 0., yeff = 0., zeff = 0., eeff = 0.;
      while(sp < (int) coreHits.size() && sp < (ihit+sparseWidth)) {
        double ecore = coreHits[sp]->Edep();
        xeff += coreHits[sp]->XYZ().X()*ecore;
        yeff += coreHits[sp]->XYZ().Y()*ecore;
        zeff += coreHits[sp]->XYZ().Z()*ecore;
        eeff += ecore;
        sp++;
      }
      ihit += sparseWidth;

      // Finish calculating effective coordinates.
      xeff /= eeff;
      yeff /= eeff;
      zeff /= eeff;

      // Create the effective hit and store it in the RBeta object.
      Hit * heff = new Hit(xeff,yeff,zeff,eeff);
      fRBeta->AddEffHit(*heff);

      echeck += eeff;
    }
    std::cout << "Total energy in sparsed hits = " << echeck << std::endl;

    // -----------------------------------------------------------------------------
    // Compute the blobs from the effective hits.
    // -----------------------------------------------------------------------------

    IBlob * blob1 = new IBlob();
    IBlob * blob2 = new IBlob();

    TVector3 cblob1 = fIBetas[0]->GetVX().first;
    TVector3 cblob2 = fIBetas[0]->GetVX().second;

    // Determine the hits in the first and second blobs given the radius.
    double eBlob1 = 0., eBlob2 = 0.;
    std::vector<const alex::Hit*> effHits = fRBeta->GetEffHits();
    for(int h = 0; h < (int) effHits.size(); h++) {

      // Store the location and energy of this hit.
      double hx = effHits[h]->XYZ().X();
      double hy = effHits[h]->XYZ().Y();
      double hz = effHits[h]->XYZ().Z();
      double he = effHits[h]->Edep();

      // Add this hit's energy to the first blob if it lies within the 
      // radius of the center.
      double dx1 = (hx - cblob1.X());
      double dy1 = (hy - cblob1.Y());
      double dz1 = (hz - cblob1.Z());
      double dist1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
      if(dist1 < blobRadius) {
        eBlob1 += he;
        blob1->AddHit(*(effHits[h]));
      }

      // Add this hit's energy to the second blob if it lies within 
      // the radius of the center.
      double dx2 = (hx - cblob2.X());
      double dy2 = (hy - cblob2.Y());
      double dz2 = (hz - cblob2.Z());
      double dist2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
      if(dist2 < blobRadius) {
        eBlob2 += he;
        blob2->AddHit(*(effHits[h]));
      }
    }

    // Add the blobs to the RBeta object.
    fRBeta->AddBlob(*blob1);
    fRBeta->AddBlob(*blob2); 

//    // Reverse the effective hits if it turned out that blob2 was actually more energetic. (THIS SHOULD ONLY BE DONE FOR "DATA")
//    if(eBlob2 > eBlob1) {
//      std::cout << log4cpp::Priority::DEBUG << " *** Reversing eff. hits because eBlob2 > eBlob1" << std::endl;
//      fRBeta->ReverseEffHits();
//    }
 
    /* NOTE: this is not trivial as the effective hits can't be easily modified after created.
    // Add all photons to their closest effective hit.
    std::vector<const IBeta*> photonTracks;
    for(int i = 0; i < (int) photonTracks.size(); i++) {
      const Hit * photCenter = photonTracks[i]->GetBlob();
      double xphot = photCenter->XYZ().X();
      double yphot = photCenter->XYZ().Y();
      double zphot = photCenter->XYZ().Z();

      // Get the effective hits.
      std::vector<const Hit*> tempHits = photonTracks[i]->GetHit();
      for(int j = 0; j < (int) tempHits.size(); j++) {
        double zhit = tempHits[j]->XYZ().Z();
        double ehit = tempHits[j]->Edep();
        if(zhit >= z && zhit <= zend) {
          eeff += ehit;
        }
      }
    }
    
    // Place all hits in a z-bin.
    std::vector<const Hit*> coreHits = fRBeta->GetCoreHits();
    for(double z = zminCore; z < zmaxCore; z += zPlaneWidth) {

      // Compute the end z-coordinate.
      double zend = z + zPlaneWidth;

      // Create an effective hit for this z.
      double xeff = 0., yeff = 0., eeff = 0.;
      for(int i = 0; i < (int) coreHits.size(); i++) {
        double zcore = coreHits[i]->XYZ().Z();
        double ecore = coreHits[i]->Edep();
        if(zcore >= z && zcore < zend) {
          xeff += coreHits[i]->XYZ().X()*ecore;
          yeff += coreHits[i]->XYZ().Y()*ecore;
          eeff += ecore;
        }
      }

      // Calculate the effective X and Y coordinates.
      xeff /= eeff;
      yeff /= eeff;

      // Add the energy of photons that are located in this z-slice.
      std::vector<const IBeta*> photonTracks;
      for(int i = 0; i < (int) photonTracks.size(); i++) {

        // Get the hits and add the energy for the hits in the corresponding z.
        std::vector<const Hit*> tempHits = photonTracks[i]->GetHit();
        for(int j = 0; j < (int) tempHits.size(); j++) {
          double zhit = tempHits[j]->XYZ().Z();
          double ehit = tempHits[j]->Edep();
          if(zhit >= z && zhit <= zend) {
            eeff += ehit;
          }
        }
      }

      // Create the effective hit and store it in the RBeta object.
      Hit * heff = new Hit(xeff,yeff,z,zend,eeff);
      fRBeta->AddEffHit(*heff);
    }*/

  }

//
// Create a vector of properly ordered effective hits and an initial
//  state (guess on position and momentum) for input to a KF.
//
// energy: the assumed energy for a single-electron track
// errXY: the x-y error
// addPhotons: true if photons within the main IBeta are to be
//  added to the hits
// forwardFit: true if a forward fit is to be performed
//--------------------------------------------------------------------
  void IreneManager::CreateKFObjects(double energy, double errXY, bool addPhotons, bool forwardFit)
//--------------------------------------------------------------------
  {

    log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::DEBUG << " IreneManager::CreateKFObjects";

    // Must have created an RBeta object before running this function.
    if(fRBeta == NULL) {
      klog << log4cpp::Priority::DEBUG << " ERROR: must call CreateRTracks before calling CreateKFObjects";
      return;
    }

    // Must have at least 2 effective hits to perform this method.
    if(fRBeta->GetEffHits().size() < 2) {
      klog << log4cpp::Priority::DEBUG << " ERROR: must have at least 2 effective hits to perform CreateKFObjects";
      return;
    }

    // Clear the current lists of hits and initial guesses.
    fKFHits.clear();
    fKFv0.clear();
    fKFp0.clear();

    // ------------------------------------------------------------------------------------
    // Create the list of properly ordered effective hits.
    // ------------------------------------------------------------------------------------

    // Create a copy of the effective hits from the rBeta object and smear them.
    std::vector<const Hit*> effHits = fRBeta->GetEffHits();
    std::vector<Hit*> updatedHits;
    for(int h = 0; h < (int) effHits.size(); h++) {

      // Smear the x and y.
      double hx = fRandom->Gaus(effHits[h]->XYZ().X(),errXY);
      double hy = fRandom->Gaus(effHits[h]->XYZ().Y(),errXY);
      double hz = effHits[h]->XYZ().Z();
      double he = effHits[h]->Edep();

      // Create and add a new hit.
      Hit * hit = new Hit(hx, hy, hz, he);
      updatedHits.push_back(hit);
    }
   
    // Add the photons to the corresponding closest hits if this option is selected.
    if(addPhotons) {
     
      std::vector<const IBeta*> photons = fRBeta->GetPhotons();
      for(int ph = 0; ph < (int) photons.size(); ph++) {

        // Iterate through the hits and determine the closest.
        const Hit * photCenter = photons[ph]->GetBlob();
        double xphot = photCenter->XYZ().X();
        double yphot = photCenter->XYZ().Y();
        double zphot = photCenter->XYZ().Z();
        double ephot = photCenter->Edep();

        // Get the effective hits and compute the distance to each.
        double mindist = 0.; int imin = -1;
        for(int j = 0; j < (int) updatedHits.size(); j++) {
          double dx = (xphot - updatedHits[j]->XYZ().X());
          double dy = (yphot - updatedHits[j]->XYZ().Y());
          double dz = (zphot - updatedHits[j]->XYZ().Z());
          double dist = sqrt(dx*dx + dy*dy + dz*dz);
          
          if(imin < 0 || dist < mindist) {
            imin = j;
            mindist = dist;
          }
        }

        // Add the photon energy to the closest hit.
        updatedHits[imin]->AddEnergy(ephot);

      }
    }

    // Save the list of hits, reversing the hits if necessary.
    if(forwardFit) {

      // Reverse the hits, so that the least energetic blob is recorded first.
      for(int h = (int) updatedHits.size()-1; h >= 0; h--)
        fKFHits.push_back(updatedHits[h]);
    }
    else {
    
      for(int h = 0; h < (int) updatedHits.size(); h++)
        fKFHits.push_back(updatedHits[h]);
    }

    // ------------------------------------------------------------------------------------
    // Compute the guesses for the initial position and momentum.
    // ------------------------------------------------------------------------------------
    
    // Use the location of the first hit as the position guess.
    double x0 = fKFHits[0]->XYZ().X();
    double y0 = fKFHits[0]->XYZ().Y();
    double z0 = fKFHits[0]->XYZ().Z();
    fKFv0.push_back(x0);
    fKFv0.push_back(y0);
    fKFv0.push_back(z0);

    klog << log4cpp::Priority::DEBUG << " Initial position guess: (" << fKFv0[0] 
         << ", " << fKFv0[1] << ", " << fKFv0[2] << ")\n";
    
    // Use the first two hits to guess the momentum direction.
    double dx0 = (fKFHits[1]->XYZ().X() - fKFHits[0]->XYZ().X());
    double dy0 = (fKFHits[1]->XYZ().Y() - fKFHits[0]->XYZ().Y());
    double dz0 = (fKFHits[1]->XYZ().Z() - fKFHits[0]->XYZ().Z());
    double mag = sqrt(dx0*dx0 + dy0*dy0 + dz0*dz0);

    double pmag = sqrt(energy*energy - 0.5109989*0.5109989);
    fKFp0.push_back(pmag*dx0/mag);
    fKFp0.push_back(pmag*dy0/mag);
    fKFp0.push_back(pmag*dz0/mag);

    /*// Perform a fit to the first few hits to guess the initial position.
    int nfit = 4;
    if(nfit > fKFHits.size()) {
      std::cout << "WARNING: not enough hits for initial direction fit with " << nfit << " points" << std::endl;
    }
    TArrayD fitData(3*nfit);
    for(int f = 0; f < nfit && f < (int) fKFHits.size(); f++) {
      fitData[3*f] = fKFHits[f]->XYZ().X() - x0;
      fitData[3*f+1] = fKFHits[f]->XYZ().Y() - y0;
      fitData[3*f+2] = fKFHits[f]->XYZ().Z() - z0;
    }
    TMatrixD fitMatrix(3,f,fitData.GetArray());*/

    klog << log4cpp::Priority::DEBUG << " Initial momentum guess: (" << fKFp0[0]
         << ", " << fKFp0[1] << ", " << fKFp0[2] << ")\n";

    // ------------------------------------------------------------------------------------
    // Set measurement errors.
    // ------------------------------------------------------------------------------------
    fKFMErrors.push_back(errXY);
    fKFMErrors.push_back(errXY);
    fKFMErrors.push_back(errXY);
    //fKFMErrors.push_back(0.);

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

    for (int i = 0; i < (int) ihs1.size(); i++)
    {
      for (int j = 0; j < (int) ihs2.size(); j++)
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
  double IreneManager::ComputeMinDist(std::vector<const alex::Hit*> ihs1, IHit ih, int& imin, int& iside)
//--------------------------------------------------------------------
  {

    // Find the minimum distance and index of the corresponding hit in the track.
    double mindist = 1.e9;
    double diffX, diffY, diffZ, dist;
    double ihx = ih.first.X();
    double ihy = ih.first.Y();
    double ihz = ih.first.Z();
    for (int i = 0; (int) i < ihs1.size(); i++)
    {
      diffX = ihs1[i]->XYZ().X() - ihx;
      diffY = ihs1[i]->XYZ().Y() - ihy;
      diffZ = ihs1[i]->XYZ().Z() - ihz;
      dist = sqrt(diffX*diffX+diffY*diffY+diffZ*diffZ);
      //std::cout << "Comparing hit (" << ih.first.X() << ", " << ih.first.Y()
      //        << ", " << ih.first.Z() << ") with hit at (X, Y, Z) = (" << ihs1[i]->XYZ().X() 
      //        << ", " << ihs1[i]->XYZ().Y() << ", " << ihs1[i]->XYZ().Z() << "): dist = " << dist
      //        << std::endl;

      if(dist < mindist) {
        mindist = dist;
        imin = i;
      }
    }

    // Determine the index (imin+1 or imin-1) closer to the hit ih.
    double distm1, distp1;
    if(ihs1.size() <= 1) iside = -1;
    else if(imin == 0) iside = imin + 1;
    else if(imin == (int) (ihs1.size()-1)) iside = imin - 1;
    else {

      // Calculate the distance to the hit at imin-1.
      diffX = ihs1[imin-1]->XYZ().X() - ihx;
      diffY = ihs1[imin-1]->XYZ().Y() - ihy;
      diffZ = ihs1[imin-1]->XYZ().Z() - ihz;
      distm1 = sqrt(diffX*diffX+diffY*diffY+diffZ*diffZ);
      
      // Calculate the distance to the hit at imin+1.
      diffX = ihs1[imin+1]->XYZ().X() - ihx;
      diffY = ihs1[imin+1]->XYZ().Y() - ihy;
      diffZ = ihs1[imin+1]->XYZ().Z() - ihz;
      distp1= sqrt(diffX*diffX+diffY*diffY+diffZ*diffZ);

      // Choose the index of the hit that is closer.
      if(distm1 < distp1) iside = imin - 1;
      else iside = imin + 1;
    }

    return mindist;
  }

  /*
//--------------------------------------------------------------------
  double IreneManager::ComputeMinDist(std::vector<const alex::Hit*> ihs1, IHit ih)
//--------------------------------------------------------------------
  {
    double mindist = 1.e9;

    for (int i = 0; i < (int) ihs1.size(); i++)
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
  */
}
