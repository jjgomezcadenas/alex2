// Generated by AlexConf: do not edit
#include "IElectrons.hh"
#include <alex/ISvc.h>
namespace alex {

//--------------------------------------------------------------------
  bool IElectrons::Init()
//--------------------------------------------------------------------
  {
  	return true;
  }
//--------------------------------------------------------------------
  bool IElectrons::Execute()
//--------------------------------------------------------------------
  {
  	log4cpp::Category& klog = log4cpp::Category::getRoot();
    klog << log4cpp::Priority::DEBUG << " IElectrons::Execute" ;


		klog << log4cpp::Priority::DEBUG << " GetElectrons()" ;
  	IParticles vBetas = ISvc::Instance().GetElectrons();
    IParticles vPrimaryBetas = ISvc::Instance().GetPrimaryElectrons();

  	auto nBetas = ISvc::Instance().GetNumberOfElectrons();
  	klog << log4cpp::Priority::DEBUG << " GetNumberElectrons() = " << nBetas;

  	fH1_nb->Fill(nBetas);

  	for (auto beta : vBetas)
  	{
  		klog << log4cpp::Priority::DEBUG << " beta->Momentum() = " << beta->Momentum();
      fH1_P->Fill(beta->Momentum());
  	}

    auto nPrimaryBetas = ISvc::Instance().GetNumberOfPrimaryElectrons();
    klog << log4cpp::Priority::DEBUG << " GetNumberPrimaryElectrons() = " 
    << nPrimaryBetas;

    fH1_npb->Fill(nPrimaryBetas);

    for (auto beta : vPrimaryBetas)
    {
      klog << log4cpp::Priority::DEBUG << " beta->Momentum() = " << beta->Momentum();
      fH1_PP->Fill(beta->Momentum());
    }

  	klog << log4cpp::Priority::DEBUG << " Calling GetPMaxElectrons()  = " ;
  	std::pair<IParticle,IParticle> bmax =ISvc::Instance().GetPMaxElectrons();

  	klog << log4cpp::Priority::DEBUG << " bmax.first->Momentum()  = "
  	<< bmax.first->Momentum();

  	fH1_PMax->Fill(bmax.first->Momentum());
  	auto TMax = bmax.first->Energy()-bmax.first->GetMass();

  	klog << log4cpp::Priority::DEBUG << " TMax  = "
  	<< TMax;

  	fH1_TMax->Fill(TMax);

  	auto TMax2 =0.;
  	if (ISvc::Instance().GetNumberOfPrimaryElectrons() >1)
  	{
  		fH1_PMax2->Fill(bmax.second->Momentum());

  		klog << log4cpp::Priority::DEBUG << " bmax.second->Momentum()  = "
  	<< bmax.second->Momentum();

  		TMax2 = bmax.second->Energy()-bmax.second->GetMass();

  		klog << log4cpp::Priority::DEBUG << " TMax2  = "
  	<< TMax2;
  		fH1_TMax2->Fill(TMax2);
  	} 

  	fH1_TMax12->Fill(TMax+TMax2);
  		

		auto vf = bmax.first->GetDecayVertex();

    klog << log4cpp::Priority::DEBUG << " decay vertex: X  = "
    <<vf[0] << " Y = " << vf[1] << " Z = " << vf[2];

		fH2_YZ->Fill(vf[1],vf[2]);
		fH2_XZ->Fill(vf[0],vf[2]);
		fH2_XY->Fill(vf[0],vf[1]);

		return true;
  }
//--------------------------------------------------------------------
  bool IElectrons::End()
//--------------------------------------------------------------------
  {
  	return true;
  }
}
