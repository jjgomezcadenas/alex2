
// ----------------------------------------------------------------------------
//  $Id: ToyAnalysis.cpp 
//
//  Author:  <gomez@mail.cern.ch>
//  Created: July 2014
//  
//  Copyright (c) 2014 NEXT Collaboration
// ---------------------------------------------------------------------------- 



#include <TCanvas.h>
#include <TStyle.h>

#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <map>

#include <alex/Alex.h>
#include <alex/ToyAnalysis.h>
#include <alex/ToyData.h>

using std::string; 
using std::cout;
using std::cin; 
using std::endl; 
using std::ostream;
using std::vector;

namespace alex {


//--------------------------------------------------------------------
  bool ToyAnalysis::Init()
//--------------------------------------------------------------------
  {
    
    fH1 = new TH1F("fH1","x distribution",100,-4,4);
    return true;

  }
//--------------------------------------------------------------------
  bool ToyAnalysis::Execute()
//--------------------------------------------------------------------
  {
    IData* itoyData =alex::Alex::Instance().RetrieveData("toyData");
    ToyData* toyData =dynamic_cast<ToyData*> (itoyData);
    TVector3 V3 = toyData->GetData();
    fH1->Fill(V3[0]);
      return true;
  }
//--------------------------------------------------------------------
  bool ToyAnalysis::End()
//--------------------------------------------------------------------
  {
    TCanvas *c1 = new TCanvas( "c1", "Hits", 200, 10, 600, 800 );
    gStyle->SetOptStat(1000000001);

    c1->cd();
    fH1->Draw();
    c1->Update();
    
    string input ;
    cout << "Return to continue:\n>";
    getline(cin, input);

    return true;
  }

}

