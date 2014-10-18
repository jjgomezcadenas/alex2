#include "BlobRadius.hh"
#include <alex/ISvc.h>
namespace alex {

  //--------------------------------------------------------------------
  bool BlobRadius::Init()
  //--------------------------------------------------------------------
  {

    return true;
  }
  //--------------------------------------------------------------------
  bool BlobRadius::Execute()
  //--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::DEBUG << " BlobRadius::Execute" ;

    IParticles vpe = ISvc::Instance().GetPrimaryElectrons();
    klog << log4cpp::Priority::DEBUG << " true electrons = " << vpe.size();
  
    klog << log4cpp::Priority::DEBUG << " fMaxDist = " << fMaxDist;
     
    //std::pair<IHits, IHits> phits = ISvc::Instance().GetPMaxElectronsHits();	
    //fH2_nth_2e->Fill(phits.first.size(), phits.second.size());
    //fH1_nth_pe->Fill(phits.first.size());
    //fH1_nth_pe->Fill(phits.second.size());
    //std::vector< IHits > tracks = CreateTracks();

    // Init the IreneManager
    ISvc::Instance().Init("DEBUG");

    // Call the ISvc function to create the separated tracks.
    ISvc::Instance().CreateTracks(fMaxDist);
    std::vector<const IBeta*> iBetas = ISvc::Instance().GetIBetas();

    // Fail if no IBeta objects were constructed.
    if(iBetas.size() == 0) {
      std::cout << "WARNING: No IBeta objects found." << std::endl;
      return false;
    }

    // Get the hits for the electrons with maximum and second maximum momentum.
    std::pair<IHits, IHits> electronhits = ISvc::Instance().GetPMaxElectronsHits();
    IHits trkFirst = electronhits.first;
    IHits trkSecond = electronhits.second;

    // Get the locations of the blob seeds (end of both tracks).
    TLorentzVector cblob1;
    TLorentzVector cblob2;
    cblob1 = trkFirst[trkFirst.size()-1].first;
    if(ISvc::Instance().GetNumberOfPrimaryElectrons() >= 2) cblob2 = trkSecond[trkSecond.size()-1].first;
    else cblob2 = trkFirst[0].first;
    
    // Fill the number of tracks histogram.
    fH1_ntrck->Fill(iBetas.size());

    // Fill the minimum distance histogram.
    std::vector<double> mdl = ISvc::Instance().GetMinDistList();
    for(int i = 0; i < (int) mdl.size(); i++)
    {
      fH1_mindist->Fill(mdl[i]);
    }

    // Find the IBeta with the greatest number of hits.
    int maxIBeta = 0, nMax = -1;
    for(int i = 0; i < (int) iBetas.size(); i++) {
      
      // Get the number of hits in this IBeta, and compare
      //  with the current maximum.
      int nHits = (int) iBetas[i]->GetHit().size();
      if(nHits > nMax) {
        maxIBeta = i;
        nMax = nHits;
      }
    }
    std::cout << "Found maximum IBeta at index " << maxIBeta << std::endl;

    // Fill the number of hits and energy histograms.
    int lenBetaMax = (int) iBetas[maxIBeta]->GetHit().size();
    double eBetaMax = iBetas[maxIBeta]->GetEnergy();
    fH1_betaNHitsMax->Fill(lenBetaMax);
    fH1_betaEMax->Fill(eBetaMax);
    for(int i = 0; i < (int) iBetas.size(); i++) {
      if(i != maxIBeta) {
        fH1_betaNHits->Fill(iBetas[i]->GetHit().size());
        fH1_betaE->Fill(iBetas[i]->GetEnergy());
      }  
    }

    // Select only events with a "long track".
    

    // Loop over the radii.
    for(double rd = fRadMin; rd < fRadMax; rd += fRadInc) {

      // Determine the energy in the first and second blob given the radius.
      double eBlob1 = 0., eBlob2 = 0.;
      std::vector<const alex::Hit*> mainTrack = iBetas[maxIBeta]->GetHit();
      for(int h = 0; h < (int) mainTrack.size(); h++) {

        // Store the location and energy of this hit.
        double hx = mainTrack[h]->XYZ().X();
        double hy = mainTrack[h]->XYZ().Y();
        double hz = mainTrack[h]->XYZ().Z();
        double he = mainTrack[h]->Edep();
 
        // Add this hit's energy to the first blob if it lies within the 
        // radius of the center.
        double dx1 = (hx - cblob1.X());
        double dy1 = (hy - cblob1.Y());
        double dz1 = (hz - cblob1.Z());
        double dist1 = sqrt(dx1*dx1 + dy1*dy1 + dz1*dz1);
        if(dist1 < rd) eBlob1 += he;

        // Add this hit's energy to the second blob if it lies within 
        // the radius of the center.
        double dx2 = (hx - cblob2.X());
        double dy2 = (hy - cblob2.Y());
        double dz2 = (hz - cblob2.Z());
        double dist2 = sqrt(dx2*dx2 + dy2*dy2 + dz2*dz2);
        if(dist2 < rd) eBlob2 += he;
      }

      // Add the results to the 2D histograms.
      fH2_radEnergy1->Fill(rd, eBlob1);
      fH2_radEnergy2->Fill(rd, eBlob2);
      fH2_blobEnergy->Fill(eBlob1,eBlob2);
    }

    return true;
  }
  //--------------------------------------------------------------------
  bool BlobRadius::End()
  //--------------------------------------------------------------------
  {
    return true;
  }

/*
  std::vector< IHits > PElectrons::CreateTracks()
  {
    std::vector< IHits > tracks;
    tracks.clear();
    IHits truehits =ISvc::Instance().GetTrueHits();

    IHit initial_ih = truehits[0];
    std::vector< IHit > initial_track;
    initial_track.push_back(initial_ih);
    truehits.erase(truehits.begin());

	

    int i = 0;
    while(truehits.size()!=0)
      {
	IHit ih = truehits[i];
	    
	double mindist = ComputeMinDist(initial_track, ih);
	if(mindist<fMaxDist)
	  {
	    initial_track.push_back(ih);
	    truehits.erase(truehits.begin());
	    i=0;
	  }
	else 
	  {	    
	    i++;	    
	  }

	if(i>=truehits.size() && truehits.size()!=0)
	  {
	    tracks.push_back(initial_track);
	    initial_track.clear();
	    initial_ih = truehits[0];
	    truehits.erase(truehits.begin());
	    initial_track.push_back(initial_ih);
	    i=0;
	  }

	    
      }


    return tracks;
  }
*/



    // std::vector< IHits > PElectrons::CreateTracks()
    // {
    //   //    log4cpp::Category& klog2 = log4cpp::Category::getRoot();
    //   //    klog2 << log4cpp::Priority::DEBUG << " CreateTracks";



    //   std::cout<<"CreateTracks"<<std::endl;

    //   std::vector< IHits > track;
    //   track.clear();
    //   std::pair<IHits, IHits> electronhits =ISvc::Instance().GetPMaxElectronsHits();
    //   IHits truehits =ISvc::Instance().GetTrueHits();

    //   //    klog2 << log4cpp::Priority::DEBUG << " CreateTracks True Hits size = "<<truehits.size();
    //   std::cout<<" CreateTracks True Hits size = "<<truehits.size()<<std::endl;


    //   ///Create Hits vector with all hits from electron 1 and 2
    //   track.push_back(electronhits.first);       
    //   for(int i = 0; i<electronhits.second.size(); i++)
    // 	{
    // 	  track[0].push_back(electronhits.second[i]);       
    // 	}


    //   //    klog2 << log4cpp::Priority::DEBUG << " CreateTracks tracks with 2 primary electrons size = "<<track.size();
    //   std::cout<<" CreateTracks tracks with 2 primary electrons size = "<<track.size()<<std::endl;



    //   for (unsigned int j=0; j<track.size(); j++)
    // 	{
    // 	  //	while(truehits.size()>=1)
    // 	  for (unsigned int l=0; l<truehits.size(); l++)

    // 	    {

    // 	      std::cout<<"track size = "<<track.size()<<" j = "<<j<<std::endl;
	    
    // 	      if(!(InVector(track[j], truehits[0])))
    // 		{
    // 		  double mindistance = ComputeMinDist(track[j], truehits[0]);
    // 		  fH1_mindist->Fill(mindistance);
    // 		  if(mindistance<fMaxDist)
    // 		    {
    // 		      track[j].push_back(truehits[0]);
    // 		      std::cout<< " CreateTracks Hit added to track. "<<truehits[0].first.X()<<", "<<mindistance<<std::endl;
    // 		      fH1_mindist_in->Fill(mindistance);
    // 		    }
    // 		  else
    // 		    {
    // 		      IHits ihs;
    // 		      std::cout<< " CreateTracks Hit added as a different track. "<<truehits[0].first.X()<<", "<<mindistance<<std::endl;

    // 		      ihs.push_back(truehits[0]);
    // 		      std::cout<< " CreateTracks Hit added to track. Track size before = "<<track.size()<<std::endl;
    // 		      track.push_back(ihs);
    // 		      std::cout<< " CreateTracks Hit added to track. Track size after = "<<track.size()<<std::endl;
    // 		      fH1_mindist_out->Fill(mindistance);
    // 		    }
    // 		}
    // 	      //	    truehits.erase(truehits.begin());
    // 	    }
    // 	}


    //   return track;

    // }



/*
    double PElectrons::ComputeMinDist(IHits ihs1, IHits ihs2)
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


    double PElectrons::ComputeMinDist(IHits ihs1, IHit ih)
    {
      double mindist = 1.e9;

      for (int i = 0; i<ihs1.size(); i++)
	{
	  double diffX = ihs1[i].first.X()-ih.first.X();
	  double diffY = ihs1[i].first.Y()-ih.first.Y();
	  double diffZ = ihs1[i].first.Z()-ih.first.Z();
	  double dist = sqrt(diffX*diffX+diffY*diffY+diffZ*diffZ);
	  if(dist<mindist)
	    mindist = dist;	  	
	}

      return mindist;
    }


    bool PElectrons::InVector(IHits ihs, IHit ih)
    {
      bool inv = false;
      for (int i = 0; i<ihs.size(); i++)
	{
	  if(ih==ihs[i])
	    inv=true;
	}
      return inv;
    }
  */
  }
