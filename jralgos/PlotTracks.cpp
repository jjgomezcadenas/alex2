#include "PlotTracks.hh"
#include <alex/ISvc.h>
namespace alex {

  //--------------------------------------------------------------------
  bool PlotTracks::Init()
  //--------------------------------------------------------------------
  {

    return true;
  }
  //--------------------------------------------------------------------
  bool PlotTracks::Execute()
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

    fEvent++;

    // Init the IreneManager
    ISvc::Instance().Init("DEBUG");

    // Call the ISvc functions to create the separated tracks and RBeta object.
    ISvc::Instance().CreateTracks(fMaxDist);
    ISvc::Instance().CreateRTracks(fCoreDist,fSparseWidth,fBlobRadius);
    ISvc::Instance().CreateKFObjects(2.468, 1.0, true, true);
    std::vector<const IBeta*> iBetas = ISvc::Instance().GetIBetas();
    RBeta * rBeta = ISvc::Instance().GetRBeta();
    std::cout << "-- Created RBeta object for event " << fEvent << " with " << iBetas.size() << " IBetas and "
              << rBeta->GetPhotons().size() << " photons " << std::endl;

    // Fail if no IBeta objects were constructed.
    if(iBetas.size() == 0) {
      std::cout << "WARNING: No IBeta objects found." << std::endl;
      return false;
    }

    if(fEvent == fEvtPlot) {

      // Plot the core hits.
      std::vector<const Hit*> coreHits = rBeta->GetCoreHits();
      for(int i = 0; i < (int) coreHits.size(); i++) {
        const Hit * h = coreHits[i];
        double x = h->XYZ().X();
        double y = h->XYZ().Y();
        fH2_xyCore->Fill(x,y);
      }

      // Plot all main-track hits.
      std::vector<const Hit*> mainTrackHits = iBetas[0]->GetHit();
      for(int i = 0; i < (int) mainTrackHits.size(); i++) {
        const Hit * h = mainTrackHits[i];
        double x = h->XYZ().X();
        double y = h->XYZ().Y();
        fH2_xyMainTrack->Fill(x,y);
      }

      // Plot the effective hits.
      std::vector<const Hit*> effHits = rBeta->GetEffHits();
      for(int i = 0; i < (int) effHits.size(); i++) {
        const Hit * h = effHits[i];
        double x = h->XYZ().X();
        double y = h->XYZ().Y();
        fH2_xyEffHits->Fill(x,y,i);
      }     

      // Plot all hits.
      IHits allHits = ISvc::Instance().GetTrueHits();
      for(int i = 0; i < (int) allHits.size(); i++) {
        double x = allHits[i].first.X();
        double y = allHits[i].first.Y();
        fH2_xyAllHits->Fill(x,y);
      }

      std::cout << "Blob1 at (" << iBetas[0]->GetVX().first.X() << ", " << iBetas[0]->GetVX().first.Y() << ", " << iBetas[0]->GetVX().first.Z() << ")" << std::endl;
      std::cout << "Blob2 at (" << iBetas[0]->GetVX().second.X() << ", " << iBetas[0]->GetVX().second.Y() << ", " << iBetas[0]->GetVX().second.Z() << ")" << std::endl;

    }  

    return true;
  }
  //--------------------------------------------------------------------
  bool PlotTracks::End()
  //--------------------------------------------------------------------
  {
    return true;
  }
}
