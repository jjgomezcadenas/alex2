
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
#include "AEx2.hh"
#include <alex/ToyData.h>

using std::string; 
using std::cout;
using std::cin; 
using std::endl; 
using std::ostream;
using std::vector;

namespace alex {


//--------------------------------------------------------------------
  bool AEx2::Init()
//--------------------------------------------------------------------
  {
  
    return true;

  }
//--------------------------------------------------------------------
  bool AEx2::Execute()
//--------------------------------------------------------------------
  {
    IData* itoyData =alex::Alex::Instance().RetrieveData("toyData");
    ToyData* toyData =dynamic_cast<ToyData*> (itoyData);
    TVector3 V3 = toyData->GetData();
    fH1_Y->Fill(V3[1]);
    fH2_YZ->Fill(V3[1],V3[2]);
      return true;
  }
//--------------------------------------------------------------------
  bool AEx2::End()
//--------------------------------------------------------------------
  {
    TCanvas *c1 = new TCanvas( "c1", "Hits1", 200, 10, 600, 800 );
    gStyle->SetOptStat(1000000001);

    c1->Divide(1,2);
    c1->cd(1);
    fH1_Y->Draw();
    c1->cd(2);
    fH2_YZ->Draw("box");
    c1->Update();
    
    string input ;
    cout << "Return to continue:\n>";
    getline(cin, input);

    return true;
  }

}

