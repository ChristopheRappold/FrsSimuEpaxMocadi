#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TGenPhaseSpace.h"
#include "TRandom3.h"
#include "TEfficiency.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "Riostream.h"
#include <vector>

void SimuPhaseSpace() 
{

  TRandom3* rand = new TRandom3;
  TH1F* h_mom_fr = new TH1F("mom_fr","mom_fr",500,6,10); 
  
  TH1F* h_InvMass = new TH1F("inv_mass","inv_mass",10000,1.5,2.5);
  TH1F* h_InvMassMix = new TH1F("inv_mass_mix","inv_mass_mix",1000,7.,8.);

  std::vector<TLorentzVector> mom_pi;
  std::vector<TLorentzVector> mom_B;
  
  int n_mix = 0;

  for(int i=0;i<100000;++i)
    {
      TLorentzVector W;
      

      //double pW = 1.9*8; // L 
      double pW = 1.8*3; // L 
      //double pW = 7.7+0.1*k; // H3L
      //double pW = 10.5+0.1*k; // H4L
      double u = 0.931494061;
      double mass  = 2.9937;//7*u+27.87e-3 + 1.11568;
  
      double pt = rand->Uniform(0.5);
      double phi = rand->Uniform(1.);
      double fractor_pW = rand->Gaus(1,0.05);
      phi *= 2*TMath::Pi();
      phi -= TMath::Pi();
      pW *= fractor_pW;
      W.SetXYZM(pt,pt,sqrt(pW*pW-pt*pt-pt*pt),mass);//,3.9225);//2.99114);
      W.SetPhi(W.Phi()+phi);

      Double_t masses[2] = {0.1396,2.80925};
      //Double_t masses[2] = {0.1396,8*u+22.9215e-3};
      
      //Double_t masses[2] = {0.1396,2.80923};
      //Double_t masses[2] = {0.1396,3.727417};
  
      TGenPhaseSpace event1,event1_2;
      
      
      //std::cout<<"Decay 1 :" <<(event1.SetDecay(W, 2, masses) ? "possible" : "not possible")<<endl;
      //std::cout<<"Decay 2 :" <<(event2.SetDecay(W, 2, masses,"Fermi") ? "possible" : "not possible")<<endl;
      
      bool status = event1.SetDecay(W, 2, masses);
      if(status==false)
	continue;

      Double_t weight = event1.Generate();
      double maxWeight = event1.GetWtMax();
      double rwt ;
      do
	{
	  weight = event1.Generate();
	  rwt = rand->Uniform(0.0, maxWeight);
	}
      while( rwt > weight );
      
      TLorentzVector *p[2];
      p[0] = event1.GetDecay(0);
      p[1] = event1.GetDecay(1);
      
      p[1]->SetXYZM(p[1]->X(),p[1]->Y(),p[1]->Z(),p[1]->M()+6.5e-3);

      const double m_d = 1.8761238;//1.875613;
      const double m_n = 0.9395653;//2.80925;
      const Double_t masses2[] = {m_n,m_d};

      status = event1_2.SetDecay(*(p[1]),2,masses2);
      //hist->h_status[i_histo]->Fill(status+2);
      if(status==false)
	return;
      
      Double_t weight1_2 = event1_2.Generate();
      Double_t maxWeight1_2 = event1_2.GetWtMax();
      do
	{
	  weight1_2 = event1_2.Generate();
	  rwt = rand->Uniform(0.0, maxWeight1_2);
	}
      while( rwt > weight1_2 );
      
      TLorentzVector *p2[2];
      p2[0] = event1_2.GetDecay(0);
      p2[1] = event1_2.GetDecay(1);
      
     

      TLorentzVector temp_p1;
      TLorentzVector temp_p2;
      double factor = rand->Gaus(1,0.02 );
      double factor2 = rand->Gaus(1,0.017);
      TVector3 temp_vec(p[0]->Vect());
      temp_p1.SetXYZM(p[0]->X()*factor,p[0]->Y()*factor,p[0]->Z()*factor,masses[0]);
      //mag = factor[j]/p[1]->P();
      temp_p2.SetXYZM(p2[1]->X()*factor2,p2[1]->Y()*factor2,p2[1]->Z()*factor2,masses2[1]);
      
      mom_pi.push_back(temp_p1);
      mom_B.push_back(temp_p2);
	  
      TLorentzVector PP;
      PP=temp_p1+temp_p2;

      if(mom_pi.size()>=n_mix)
	{
	  for(unsigned int ii=0;ii<mom_pi.size();++ii)
	    for(unsigned int jj=0;jj<mom_B.size();++jj)
	      if(ii!=jj)
		{

		  TLorentzVector mix_PP;
		  mix_PP=mom_pi[ii]+mom_B[ii];
		  h_InvMassMix->Fill(mix_PP.M());
		}
	  
	  mom_pi.clear();
	  mom_B.clear();
	}

      //double angle = p[0]->Angle(p[1]->Vect());
      //h2[k]->Fill(angle/TMath::Pi()*180.,p[0]->P(),weight*weight_eff[0]);
      //h1[k]->Fill(angle/TMath::Pi()*180.,p[0]->P());
      
      //TLorentzVector MM;
      //MM=*p[0]+*p[1];
      //h_c->Fill(k,MM.M(),weight_eff[0]*weight);

      h_mom_fr->Fill(p[1]->Pz());
      h_InvMass->Fill(PP.M());
    }

  TCanvas* c1 = new TCanvas("c1","c1",500,500);
  c1->cd();
  h_mom_fr->Draw();
  TCanvas* c2 = new TCanvas("c2","c2",500,500);
  //c2->Divide(2);
  c2->cd(1);
  h_InvMass->Draw();
  // c2->cd(2);
  // h_InvMass->Draw();
  // h_InvMassMix->Draw("same");


}
