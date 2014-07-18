
// ----------------------------------------------------------------------------
//  $Id: ToyAnalysis.cpp 
//
//  Author:  <gomez@mail.cern.ch>
//  Created: July 2014
//  
//  Copyright (c) 2014 NEXT Collaboration
// ---------------------------------------------------------------------------- 


#include <Riostream.h>
#include <TCanvas.h>
#include <TStyle.h>

#include <string>
#include <vector>
#include <utility>
#include <memory>
#include <map>
#include <alex/LogUtil.h>
#include <alex/ToyAnalysis.h>

using std::string; 
using std::cout; 
using std::endl; 
using std::ostream;
using std::vector;

namespace alex {


//--------------------------------------------------------------------
  bool  ToyAnalysis::AlexManager::Init()
//--------------------------------------------------------------------
  {
    fName="ToyAnalysis";
    
    fIn.open(Form("%sbasic.dat");
    fH1 = new TH1F("fH1","x distribution",100,-4,4);

  }
//--------------------------------------------------------------------
  void ToyAnalysis::AlexManager::Execute()
//--------------------------------------------------------------------
  {
    while (1) {
      fIn >> x >> y >> z;
      if (!fIn.good()) break;
      if (nlines < 5) printf("x=%8f, y=%8f, z=%8f\n",x,y,z);
      fH1->Fill(x);
      nlines++;
   }
  }
//--------------------------------------------------------------------
  void ToyAnalysis::AlexManager::End()
//--------------------------------------------------------------------
  {
    TCanvas *c1 = new TCanvas( "c1", "Hits", 200, 10, 600, 800 );
    gStyle->SetOptStat(1000000001);

    c1->cd();
    fH1->Draw();
    
    string input ;
    cout << "Return to continue:\n>";
    getline(cin, input);
  }

}

