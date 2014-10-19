#include "NDrawEffHits.hh"
#include <alex/ISvc.h>
#include <TH2F.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TApplication.h>
namespace alex {

  //--------------------------------------------------------------------
  bool NDrawEffHits::Init()
  //--------------------------------------------------------------------
  {
    // Init the IreneManager
    ISvc::Instance().Init("DEBUG");
  
    char * ch = "";
    int argc = 0;
    TApplication * gApplication = new TApplication("NDrawEffHits", &argc, &ch);
    gApplication->SetReturnFromRun(true);

    //log4cpp::Category& klog = log4cpp::Category::getRoot(); 
    //app->Run(kTRUE);

    return true;
  }
  //--------------------------------------------------------------------
  bool NDrawEffHits::Execute()
  //--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::DEBUG << " NDrawEffHits::Execute";

    // Construct histograms containing the effective hits and all hits
    //  if the draw option is on.
    RBeta * rBeta = ISvc::Instance().GetRBeta(); 

    if(fDraw == 1) {

      char userInput;
      TCanvas * c1 = new TCanvas("c1", "Hits", 10, 10, 600, 600);
      TH2F* h2_xyCore = new TH2F("h2_xyCore","Effective hits",200,-100.,100.,200,-100.,100.);
      TH2F* h2_xyAll = new TH2F("h2_xyAll","All hits",200,-100.,100.,200,-100.,100.);

      // Plot the core hits.
      std::vector<const Hit*> effHits = rBeta->GetEffHits();
      for(int i = 0; i < (int) effHits.size(); i++) {
        const Hit * h = effHits[i];
        double x = h->XYZ().X();
        double y = h->XYZ().Y();
        h2_xyCore->Fill(x,y,i);
      }

      // Plot all hits.
      // IHits allHits = ISvc::Instance().GetTrueHits();
      // for(int i = 0; i < (int) allHits.size(); i++) {
      //   double x = allHits[i].first.X();
      //   double y = allHits[i].first.Y();
      //   h2_xyAll->Fill(x,y);
      // }

      // Draw the histogram.
      std::cout << " Drawing histogram..." << std::endl;
      h2_xyCore->GetXaxis()->SetTitle("x");
      h2_xyCore->GetYaxis()->SetTitle("y");
      h2_xyCore->Draw("colz");
      // h2_xyAll->SetMarkerStyle(6);
      // h2_xyAll->Draw("SAME");
      c1->Update();

      // Process ROOT events allowing for the drawing of the histogram.
      //gSystem->ProcessEvents();

      // Await input from the user before continuing.
      std::cin >> userInput;

      // Close the plot and delete the histogram.
      c1->Close();
      delete h2_xyCore;
      delete h2_xyAll;
      delete c1;
    }

    return true;
  }
  //--------------------------------------------------------------------
  bool NDrawEffHits::End()
  //--------------------------------------------------------------------
  {
    return true;
  }
}
