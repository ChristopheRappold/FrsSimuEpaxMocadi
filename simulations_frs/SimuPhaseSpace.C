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
#include "THStack.h"
#include "TMultiGraph.h"

#include <vector>
#include <map>
#include <string>

void SimuPhaseSpace(int type,int TwoOrThree = 0, double rapidity_center=0.0, double scaling =0.02) 
{
  TRandom3* rand = new TRandom3;

  std::map<int,double> MapMassMother;
  MapMassMother.insert(std::pair<int,double>(0,2.9937));
  MapMassMother.insert(std::pair<int,double>(1,2.99114));
  MapMassMother.insert(std::pair<int,double>(2,3.9225));
  MapMassMother.insert(std::pair<int,double>(3,11.35871));
  MapMassMother.insert(std::pair<int,double>(4,1.116));
  MapMassMother.insert(std::pair<int,double>(5,1.19744));
  MapMassMother.insert(std::pair<int,double>(6,1.3217));

  std::map<int,double> MapChargeDaugthers;
  MapChargeDaugthers.insert(std::pair<int,double>(0,1.));
  MapChargeDaugthers.insert(std::pair<int,double>(1,2.));
  MapChargeDaugthers.insert(std::pair<int,double>(2,2.));
  MapChargeDaugthers.insert(std::pair<int,double>(3,7.));
  MapChargeDaugthers.insert(std::pair<int,double>(4,1.));
  MapChargeDaugthers.insert(std::pair<int,double>(5,1.));
  MapChargeDaugthers.insert(std::pair<int,double>(6,1.));
  
  
  std::map<int,std::vector<double> > MapMassDaugthers;
  std::vector<double> temp_vec(2);
  
  temp_vec[0] = 0.1396;
  temp_vec[1] = 2.80925;
  MapMassDaugthers.insert(std::pair<int,std::vector<double> >(0,temp_vec));

  temp_vec[0] = 0.1396;
  temp_vec[1] = 2.80923;
  MapMassDaugthers.insert(std::pair<int,std::vector<double> >(1,temp_vec));

  temp_vec[0] = 0.1396;
  temp_vec[1] = 3.727417;
  MapMassDaugthers.insert(std::pair<int,std::vector<double> >(2,temp_vec));
			  
  temp_vec[0] = 0.1396;
  temp_vec[1] = 11.1917;
  MapMassDaugthers.insert(std::pair<int,std::vector<double> >(3,temp_vec));

  temp_vec[0] = 0.1396;
  temp_vec[1] = 0.938;
  MapMassDaugthers.insert(std::pair<int,std::vector<double> >(4,temp_vec));

  temp_vec[0] = 0.1396;
  temp_vec[1] = 0.939565;
  MapMassDaugthers.insert(std::pair<int,std::vector<double> >(5,temp_vec));

  temp_vec[0] = 0.1396;
  temp_vec[1] = 1.1157;
  MapMassDaugthers.insert(std::pair<int,std::vector<double> >(6,temp_vec));

  std::map<int,std::vector<double> > MapInvMass;
  temp_vec[0] = TwoOrThree == 1 ? 1.5 : 2.7;
  temp_vec[1] = TwoOrThree == 1 ? 2.5 : 3.2;
  MapInvMass.insert(std::pair<int,std::vector<double> >(0,temp_vec));
  temp_vec[0] = 2.7;
  temp_vec[1] = 3.2;
  MapInvMass.insert(std::pair<int,std::vector<double> >(1,temp_vec));
  temp_vec[0] = 3.7;
  temp_vec[1] = 4.2;
  MapInvMass.insert(std::pair<int,std::vector<double> >(2,temp_vec));
  temp_vec[0] = 11.2;
  temp_vec[1] = 11.7;
  MapInvMass.insert(std::pair<int,std::vector<double> >(3,temp_vec));
  temp_vec[0] = 1.0;
  temp_vec[1] = 1.5;
  MapInvMass.insert(std::pair<int,std::vector<double> >(4,temp_vec));
  temp_vec[0] = 1.0;
  temp_vec[1] = 1.5;
  MapInvMass.insert(std::pair<int,std::vector<double> >(5,temp_vec));
  temp_vec[0] = 1.2;
  temp_vec[1] = 1.7;
  MapInvMass.insert(std::pair<int,std::vector<double> >(6,temp_vec));

  double KinematicsMother[7][2] = { {5.833,0.2323},{8.497,0.3138},{11.155,0.33704},{33.3526,1.00},{5.51,0.0000020},{/*6.62*/0.,0.00000020},{7.31,0.00000020}};
  double Brho_range [7][2] = { {12.,22.},{7.,17},{12.,22.},{10,20},{2,20},{2,20},{2,20}};
  if(TwoOrThree==0)
    {
      Brho_range[0][0] = 20;
      Brho_range[0][1] = 30;
    }
  if(rapidity_center>0.1)
    {
      Brho_range[0][0] = 10;
      Brho_range[0][1] = 30;
    }
  double Invmass_MinMax [] = {MapInvMass[type][0],MapInvMass[type][1]};

  TH1F* h_mom_fr = new TH1F("mom_fr","mom_fr",500,0,5); 
  
  TH1F* h_InvMass = new TH1F("inv_mass","inv_mass",10000,Invmass_MinMax[0],Invmass_MinMax[1]);
  TH1F* h_InvMassMix = new TH1F("inv_mass_mix","inv_mass_mix",10000,Invmass_MinMax[0],Invmass_MinMax[1]);

  TH2F* h_momFr_InvMass = new TH2F("momFrInv","momFrInv",1000,Brho_range[type][0],Brho_range[type][1],10000,Invmass_MinMax[0],Invmass_MinMax[1]);
  TH2F* h_momFr_momPi = new TH2F("momFr_momPi","momFr_momPi",1000,Brho_range[type][0],Brho_range[type][1],500.,0,1.5);
  TH2F* h_BrhoFr_BrhoPi = new TH2F("BrhoFr_BrhoPi","BrhoFr_BrhoPi",1000,Brho_range[type][0],Brho_range[type][1],500.,0,5);
  
  TH2F* h_PzFr_InvMass = new TH2F("PzFrInv","PzFrInv",1000,Brho_range[type][0],Brho_range[type][1],10000,Invmass_MinMax[0],Invmass_MinMax[1]);
  TH2F* h_PzFr_InvMassMix = new TH2F("PzFrInvMix","PzFrInvMix",1000,Brho_range[type][0],Brho_range[type][1],10000,Invmass_MinMax[0],Invmass_MinMax[1]);

  TH2F* h_MomThetaPi = new TH2F("MomThetaPi","MomThetaPi",360,0,180,300,0,3);
  TH2F* h_MomThetaFr = new TH2F("MomThetaFr","MomThetaFr",360,0,180,1800,0,3*18);
  
  THStack* h_all = new THStack("h_all","h_all");
  THStack* h_all2D = new THStack("h_all2D","h_all2D");

  TH2F* h_phasespace = new TH2F("h_phasespace","h_phasespace",100,0,4,100,0,1);
  
  std::vector<TLorentzVector> mom_pi;
  std::vector<TLorentzVector> mom_B;
  
  unsigned int n_mix = 4;

  double mass_mother = MapMassMother[type];
  double charge_D2 = MapChargeDaugthers[type];


  for(int i=0;i<100000;++i)
    {
      TLorentzVector W;
      
      // double beta_W = rand->Uniform(0.8,0.99);
      // double gamma_W = 1./TMath::Abs(1-beta_W*beta_W)
      
      //double pW = 1.9*8; // L 
      double pW = rand->Gaus(KinematicsMother[type][0],KinematicsMother[type][1]);// rand->Gaus(5.833,0.2323);//rand->Uniform(1.8,2.1)*3; // L 
      //double pW = 7.7+0.1*k; // H3L
      //double pW = 10.5+0.1*k; // H4L
      double u = 0.931494061;
      double mass  = mass_mother;//2.9937;//7*u+27.87e-3 + 1.11568;
  
      double pt = rand->Uniform(0.001);
      double phi = rand->Uniform(1.);
      double fractor_pW = 1.0;//rand->Gaus(1,0.05);
      phi *= 2*TMath::Pi();
      phi -= TMath::Pi();
      pW *= fractor_pW;
      W.SetXYZM(pt,pt,pW,mass);//sqrt(pW*pW-pt*pt-pt*pt),mass);//,3.9225);//2.99114);
      W.SetPhi(W.Phi()+phi);
      //double y = W.Rapidity();
      // double beta_W = TMath::tanH(y);
      if(rapidity_center>0.01 || type == 0)
	{
	  if(type==0)
	    {
	      double beta_Diff;
	      if(rapidity_center>0.01)
		 beta_Diff = TMath::TanH(rapidity_center-1.45);
	      else
		beta_Diff = TMath::TanH(1.8-1.45);
	      W.Boost(0,0,beta_Diff);
	    }
	  else
	    {
	      double beta_Diff = TMath::TanH(rapidity_center-1.8);
	      W.Boost(0,0,beta_Diff);
	    }
	}
      Double_t masses[2] = {MapMassDaugthers[type][0],MapMassDaugthers[type][1]};//{0.1396,2.80925};
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
      
      TLorentzVector *p2[2];
      double mass_d2 = 0;
      if(TwoOrThree==1)
	{
	  p[1]->SetXYZM(p[1]->X(),p[1]->Y(),p[1]->Z(),p[1]->M()+6.5e-3);
	  
	  const double m_d = 1.8761238;//1.875613;
	  const double m_n = 0.9395653;//2.80925;
	  const Double_t masses2[] = {m_n,m_d};
	  mass_d2 = m_d;
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
	  
	  p2[0] = event1_2.GetDecay(0);
	  p2[1] = event1_2.GetDecay(1);
	}
      else
	{
	  p2[1] = p[1];// = event1.GetDecay(1);
	  mass_d2 = masses[1];
	}
     

      TLorentzVector temp_p1;
      TLorentzVector temp_p2;
      double factor = rand->Gaus(1,1.e-2 );//0.02);
      double factor2 = rand->Gaus(1,1.e-3);//0.017);
      TVector3 temp_vec(p[0]->Vect());
      temp_p1.SetXYZM(p[0]->X()*factor,p[0]->Y()*factor,p[0]->Z()*factor,masses[0]);
      //mag = factor[j]/p[1]->P();
      temp_p2.SetXYZM(p2[1]->X()*factor2,p2[1]->Y()*factor2,p2[1]->Z()*factor2,mass_d2);
      
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
		  mix_PP=mom_pi[ii]+mom_B[jj];
		  h_InvMassMix->Fill(mix_PP.M());
		  h_PzFr_InvMassMix->Fill(3.10715497*mom_B[jj].P()/charge_D2,mix_PP.M());
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

      h_mom_fr->Fill((W.E()-W.M())/(double)(TMath::Nint(mass_mother)));//p[1]->Pz());
      h_phasespace->Fill(W.Rapidity(),W.Pt());
      h_InvMass->Fill(PP.M());
      h_momFr_momPi->Fill(3.30715497*p2[1]->Pz()/charge_D2,p[0]->P());
      h_BrhoFr_BrhoPi->Fill(3.30715497*p2[1]->Pz()/charge_D2,3.30715497*p[0]->P());
      h_momFr_InvMass->Fill(3.30715497*p2[1]->P()/charge_D2,PP.M());
      h_PzFr_InvMass->Fill(3.30715497*p2[1]->P()/charge_D2,PP.M());
      h_MomThetaFr->Fill(p2[1]->Theta()*TMath::RadToDeg(),p2[1]->P());
      h_MomThetaPi->Fill(p[0]->Theta()*TMath::RadToDeg(),p[0]->P());
    }


  
  TH2F* h_3 = (TH2F*)h_PzFr_InvMass->Clone("h_BrhoInv2");
  double binX_3 = h_3->GetXaxis()->GetBinWidth(1);
  int rebin3 = 0.008/binX_3;
  cout<<"bin X "<<binX_3<<" "<<rebin3<<" "<<0.008/binX_3<<endl;

  //h_3->Rebin2D(4,1);
  
  TH2F* h_4 = (TH2F*)h_PzFr_InvMassMix->Clone("h_BrhoInvMix2");
  double binX_4 = h_4->GetXaxis()->GetBinWidth(1);
  //int rebin4 = 0.008/binX_4;
  //h_4->Rebin2D(4,1);
  
  TGraph* a1_BrhoInvPeak = new TGraph;
  TGraph* a1_BrhoInvMin = new TGraph;
  TGraph* a1_BrhoInvMax = new TGraph;
  TGraph* a1_BrhoInvMixPeak = new TGraph;
  TGraph* a1_BrhoInvMixMin = new TGraph;
  TGraph* a1_BrhoInvMixMax = new TGraph;

  a1_BrhoInvPeak->SetNameTitle("InvPeak","InvPeak"); 
  a1_BrhoInvMin->SetNameTitle("InvMin","InvMin"); 
  a1_BrhoInvMax->SetNameTitle("InvMax","InvMax"); 
  a1_BrhoInvMixPeak->SetNameTitle("InvMixPeak","InvMixPeak"); 
  a1_BrhoInvMixMin->SetNameTitle("InvMixMin","InvMixMin"); 
  a1_BrhoInvMixMax->SetNameTitle("InvMixMax","InvMixMax"); 

  TH2F* a1_BrhoInvRMS = new TH2F("Brho_InvRMS","Brho_InvRMS",1000,Brho_range[type][0],Brho_range[type][1],1000,0.01,0.001);
  TH2F* a1_BrhoInvMixRMS = new TH2F("Brho_MixRMS","Brho_MixRMS",1000,Brho_range[type][0],Brho_range[type][1],1000,0.01,0.001);

  TH2F* a1_BrhoInvM = new TH2F("Brho_InvMean","Brho_InvMean",1000,Brho_range[type][0],Brho_range[type][1],500,Invmass_MinMax[0],Invmass_MinMax[1]);
  TH2F* a1_BrhoInvMixM = new TH2F("Brho_MixMean","Brho_MixMean",1000,Brho_range[type][0],Brho_range[type][1],500,Invmass_MinMax[0],Invmass_MinMax[1]);
  
  for(int ii=1;ii<=h_3->GetXaxis()->GetNbins();++ii)
    {
      TString nametemp("hist_proj");
      nametemp+=ii;
      TH1D* h_temp = h_3->ProjectionY(nametemp,ii,ii);
      double brho_temp = h_3->GetXaxis()->GetBinCenter(ii);
      // cout<<nametemp<<" "<<brho_temp<<" "<<h_temp->GetEntries();
      if(h_temp->GetEntries()>10)
	{
	  double mean = h_temp->GetMean();
	  double sigma = h_temp->GetRMS();
	  int binMax = h_temp->GetMaximumBin();
	  double max = h_temp->GetBinCenter(binMax);
	  //cout<<" "<<mean<<" "<<sigma;
	  //a1_BrhoInv->SetPoint(a1_BrhoInv->GetN(),brho_temp,sigma);
	  a1_BrhoInvRMS->Fill(brho_temp,sigma);
	  a1_BrhoInvM->Fill(brho_temp,max);
	  if(ii%15==0)
	    {
	      a1_BrhoInvPeak->SetPoint(a1_BrhoInvPeak->GetN(),brho_temp,max);
	      a1_BrhoInvMin->SetPoint(a1_BrhoInvMin->GetN(),brho_temp,max-3.*sigma);
	      a1_BrhoInvMax->SetPoint(a1_BrhoInvMax->GetN(),brho_temp,max+3.*sigma);
	    }
	}
      h_temp->Delete();
      // cout<<endl;
    }

  for(int ii=1;ii<=h_4->GetXaxis()->GetNbins();++ii)
    {
      TString nametemp("histMix_proj");
      nametemp+=ii;
      TH1D* h_temp = h_4->ProjectionY(nametemp,ii,ii);
      double brho_temp = h_4->GetXaxis()->GetBinCenter(ii);
      if(h_temp->GetEntries()>10)
	{
	  double mean = h_temp->GetMean();
	  double sigma = h_temp->GetRMS();
	  int binMax = h_temp->GetMaximumBin();
	  double max = h_temp->GetBinCenter(binMax);
	  //a1_BrhoInvMix->SetPoint(a1_BrhoInvMix->GetN(),brho_temp,sigma);
	  a1_BrhoInvMixRMS->Fill(brho_temp,sigma);
	  a1_BrhoInvMixM->Fill(brho_temp,max);
	  if(ii%15==0)
	    {
	      a1_BrhoInvMixPeak->SetPoint(a1_BrhoInvMixPeak->GetN(),brho_temp,mean);
	      a1_BrhoInvMixMin->SetPoint(a1_BrhoInvMixMin->GetN(),brho_temp,mean-3.*sigma);
	      a1_BrhoInvMixMax->SetPoint(a1_BrhoInvMixMax->GetN(),brho_temp,mean+3.*sigma);
	    }
	}
      h_temp->Delete();
    }


  TCanvas* c1 = new TCanvas("c1","c1",1000,500);
  c1->Divide(2);
  c1->cd(1);
  h_mom_fr->Draw();
  c1->cd(2);
  h_phasespace->Draw("colz");
  TCanvas* c2 = new TCanvas("c2","c2",1000,500);
  c2->Divide(2);
  c2->cd(1);
  h_InvMass->Draw();
  c2->cd(2);
  // h_InvMass->Draw();
  h_InvMassMix->Draw();
  TCanvas* c3_0 = new TCanvas("c3_0","c3_0",500,500);
  //c3->Divide(2,2);
  c3_0->cd(1);
  h_PzFr_InvMass->Draw("colz");
  
  TCanvas* c3_1 = new TCanvas("c3_1","c3_1",500,500);
  c3_1->cd(1);
  h_PzFr_InvMassMix->Draw("colz");

  TCanvas* c3_2 = new TCanvas("c3_2","c3_2",500,500);
  c3_2->cd(1);
  h_momFr_momPi->Draw("colz");

  TCanvas* c3_22 = new TCanvas("c3_22","c3_22",500,500);
  c3_22->cd(1);
  h_BrhoFr_BrhoPi->Draw("colz");

  
  TCanvas* c3_3 = new TCanvas("c3_3","c3_3",500,500);
  c3_3->cd(1);
  TH2F* hh_1 = (TH2F*)  h_PzFr_InvMass->Clone("brhoInv");
  TH2F* hh_2 = (TH2F*)  h_PzFr_InvMassMix->Clone("brhoMix");
  hh_1->Rebin2D(40,50);
  hh_2->Rebin2D(40,50);
  hh_2->SetLineColor(1);
  hh_2->Draw("candle");
  hh_1->SetLineColor(2);
  hh_1->Draw("box same");

  TCanvas* c4 = new TCanvas("c4","c4",500,500);
  c4->cd();
  TH1F* h_1 = (TH1F*)h_InvMass->Clone("h_Inv2");
  h_1->RebinX(10);
  h_1->Scale(scaling);
  double mean_1 = h_1->GetMean();
  double rms_1 = h_1->GetRMS();
  h_1->GetXaxis()->SetRangeUser(mean_1-3.*rms_1,mean_1+3.*rms_1);
  double integral1 = h_1->Integral();
  h_1->GetXaxis()->SetRangeUser(mean_1-30.*rms_1,mean_1+30.*rms_1);
  h_1->SetLineColor(kBlack);
  h_1->SetFillColor(kRed);
  h_1->SetLineWidth(1);
  TH1F* h_2 = (TH1F*)h_InvMassMix->Clone("h_Mix2");
  h_2->SetLineColor(kBlack);
  h_2->SetLineWidth(1);
  h_2->RebinX(10);
  h_2->GetXaxis()->SetRangeUser(mean_1-3.*rms_1,mean_1+3.*rms_1);
  double integral2 = h_2->Integral();
  h_2->GetXaxis()->SetRangeUser(mean_1-300.*rms_1,mean_1+300.*rms_1);
  cout<<" S/B :"<<integral1/integral2;
  h_1->SetLineWidth(3);
  h_2->SetLineWidth(3);
  h_all->Add(h_2);
  h_all->Add(h_1);
  h_all->Draw("");


  TCanvas* c4_1 = new TCanvas("c4_1","c4_1",500,500);
  c4_1->cd();
  TH2F* hh_4 = (TH2F*)h_PzFr_InvMass->Clone("h_BrInv2");
  hh_4->Rebin2D(5,10);
  hh_4->Scale(0.1);

  hh_4->SetLineColor(kBlack);
  hh_4->SetFillColor(kRed);
  hh_4->SetLineWidth(1);
  TH2F* hh_5 = (TH2F*)h_PzFr_InvMassMix->Clone("h_BrMix2");
  hh_5->SetLineColor(kBlack);
  hh_5->SetLineWidth(1);
  hh_5->Rebin2D(5,10);
  
  h_all2D->Add(hh_4);
  h_all2D->Add(hh_5);
  h_all2D->Draw("box");


  TCanvas* c5_1 = new TCanvas("c5_1","c5_1",500,500);
  //c5->Divide(2,2);
  c5_1->cd(1);
  //a1_BrhoInv->Draw("a*");
  a1_BrhoInvRMS->Rebin2D(10,1);
  a1_BrhoInvRMS->Draw("box");

  TCanvas* c5_2 = new TCanvas("c5_2","c5_2",500,500);
  c5_2->cd(1);
  //a1_BrhoInvMix->Draw("a*");
  a1_BrhoInvMixRMS->Rebin2D(10,1);
  a1_BrhoInvMixRMS->Draw("box");
  
  TCanvas* c5_3 = new TCanvas("c5_3","c5_3",500,500);
  c5_3->cd(1);
  a1_BrhoInvM->Draw("box");
  
  TCanvas* c5_4 = new TCanvas("c5_4","c5_4",500,500);
  c5_4->cd(1);
  a1_BrhoInvMixM->Rebin2D(20,2);
  a1_BrhoInvMixM->SetLineColor(kBlack);
  a1_BrhoInvM->Rebin2D(20,2);
  a1_BrhoInvM->SetLineColor(kRed);
  a1_BrhoInvMixM->Draw("box");
  a1_BrhoInvM->Draw("box same");

  TCanvas* c5_5 = new TCanvas("c5_5","c5_5",500,500);
  c5_5->cd(1);
  TMultiGraph* mg = new TMultiGraph();
  a1_BrhoInvMixPeak->SetMarkerColor(kRed); 
  a1_BrhoInvMixMin->SetMarkerColor(kRed); 
  a1_BrhoInvMixMin->SetLineColor(kRed); 
  a1_BrhoInvMixMax->SetMarkerColor(kRed);
  a1_BrhoInvMixMax->SetLineColor(kRed); 
  a1_BrhoInvMixPeak->SetMarkerStyle(4);
  a1_BrhoInvMixMin->SetMarkerStyle(6); 
  a1_BrhoInvMixMax->SetMarkerStyle(6); 
 
  mg->Add(a1_BrhoInvMixPeak,"p");
  mg->Add(a1_BrhoInvMixMin,"lp");
  mg->Add(a1_BrhoInvMixMax,"lp");
  // a1_BrhoInvMixMin->Draw("a*");
  // a1_BrhoInvMixMax->Draw("*");

  // a1_BrhoInvMin->Draw("*");
  // a1_BrhoInvMax->Draw("*");
  a1_BrhoInvPeak->SetMarkerColor(kBlue); 
  a1_BrhoInvMin->SetMarkerColor(kBlue); 
  a1_BrhoInvMin->SetLineColor(kBlue); 
  a1_BrhoInvMax->SetMarkerColor(kBlue);
  a1_BrhoInvMax->SetLineColor(kBlue);
  a1_BrhoInvPeak->SetMarkerStyle(8); 
  a1_BrhoInvMin->SetMarkerStyle(6); 
  a1_BrhoInvMax->SetMarkerStyle(6); 

  
  mg->Add(a1_BrhoInvPeak,"p");
  mg->Add(a1_BrhoInvMin,"lp");
  mg->Add(a1_BrhoInvMax,"lp");

  mg->Draw("a");

  TCanvas* c6 = new TCanvas("c6","c6",1000,500);
  c6->Divide(2,1);
  c6->cd(1);
  h_MomThetaPi->Draw("colz");
  c6->cd(2);
  h_MomThetaFr->Draw("colz");
  
  
  
}
