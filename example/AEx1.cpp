
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
#include "AEx1.hh"
#include "ExData.h"
#include <alex/LogUtil.h>

using std::string; 
using std::cout;
using std::cin; 
using std::endl; 
using std::ostream;
using std::vector;

namespace alex {


//--------------------------------------------------------------------
  bool AEx1::Init()
//--------------------------------------------------------------------
  {
    SetDebugLevel("INFO");
    return true;

  }
//--------------------------------------------------------------------
  bool AEx1::Execute()
//--------------------------------------------------------------------
  {
    log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::DEBUG 
        << " in AEx1::Execute, now Retrieve Data " ;

    IData* itoyData =alex::Alex::Instance().RetrieveData("ExData");

    if (itoyData == NULL)
    {
      klog << log4cpp::Priority::ERROR << "Error retrieveing data";
      exit(-1);
    }
    ExData* toyData =dynamic_cast<ExData*> (itoyData);
    TVector3 V3 = toyData->GetData();

    klog << log4cpp::Priority::DEBUG 
        << " V3[0] = " << V3[0]
        << " V3[1] = " << V3[1]
        << " V3[2] = " << V3[2];

    fH1_X->Fill(V3[0]);
    fH2_XY->Fill(V3[0],V3[1]);
    return true;
  }
//--------------------------------------------------------------------
  bool AEx1::End()
//--------------------------------------------------------------------
  {
    TCanvas *c1 = new TCanvas( "c1", "Hits", 200, 10, 600, 800 );
    gStyle->SetOptStat(1000000001);

    c1->Divide(1,2);
    c1->cd(1);
    fH1_X->Draw();
    c1->cd(2);
    fH2_XY->Draw("box");
    c1->Update();
    
    string input ;
    cout << "Return to continue:\n>";
    getline(cin, input);

    return true;
  }

}

