#include "Riostream.h"
#include "TGraph.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include <sstream>

#include <map>

#include "TRandom3.h"

#include "TH1F.h"
#include "TF1.h"

#ifdef UNFOLD
#include "../RooUnfold-1.1.1/src/RooUnfoldResponse.h"
#include "../RooUnfold-1.1.1/src/RooUnfoldBayes.h"
#include "../RooUnfold-1.1.1/src/RooUnfoldSvd.h"
#include "../RooUnfold-1.1.1/src/RooUnfoldTUnfold.h"
#endif

void LoadResEff(std::string namefile)
{

  std::ifstream ifs ( namefile.c_str() );
  std::string temp_line;
  std::vector<std::map<double,double > > energy_eff(3);

  std::getline(ifs,temp_line);

  while(std::getline(ifs,temp_line))
    {
      std::stringstream stream(temp_line);
      double energy,S,eff,eff_tot;
     
      stream >> energy >> S >> eff >> eff_tot;
      //std::cout<<" S:"<<S<<" "<<energy<<" "<<eff<<" "<<eff_tot<<std::endl;
      energy_eff[S-2].insert(std::pair<double,double>(energy,eff));

    }

  TF1* f_s = new TF1("f1","[0]*TMath::Gaus(x,[1],[2],1)",.3,3.4);
  f_s->SetParameter(0,1);
  f_s->SetParameter(1,1.805088);
  f_s->SetParameter(2,0.1801583);

  TGraph* a1_S3 = new TGraph;
  a1_S3->SetNameTitle("S3_eff","S3_eff");
  TGraph* a1_S4 = new TGraph;
  a1_S4->SetNameTitle("S4_eff","S4_eff");

  for(std::map<double,double>::iterator it_energy = energy_eff[1].begin(),it_energy_end = energy_eff[1].end();it_energy!=it_energy_end;++it_energy)
    {
      std::map<double,double>::iterator it_eff_S2 = energy_eff[0].find(it_energy->first);
      double eff_S2 = it_eff_S2->second;
      double eff_S3 = it_energy->second;

      double int_f_s = f_s->Integral(it_energy->first*1.e-3-0.01,it_energy->first*1.e-3+0.01);
      
      a1_S3->SetPoint(a1_S3->GetN(),it_energy->first*1.e-3,eff_S3/eff_S2/int_f_s);
    }

  for(std::map<double,double>::iterator it_energy = energy_eff[2].begin(),it_energy_end = energy_eff[2].end();it_energy!=it_energy_end;++it_energy)
    {
      std::map<double,double>::iterator it_eff_S2 = energy_eff[0].find(it_energy->first);
      double eff_S2 = it_eff_S2->second;
      double eff_S4 = it_energy->second;

      double int_f_s = f_s->Integral(it_energy->first*1.e-3-0.01,it_energy->first*1.e-3+0.01);
      
      a1_S4->SetPoint(a1_S4->GetN(),it_energy->first*1.e-3,eff_S4/eff_S2/int_f_s);
    }

  //TF1* f_s = new TF1("f1","TMath::Gaus(x,[0],[1],[2])",.3,3.4);
  
  int nDim = a1_S3->GetN();


  TH1F* h1_s3 = new TH1F("h1_s3","h1_s3",50,1.414,2.414);
  //double max_f = f_s->GetMaximum(1.4,2.4);
  for(int i=0;i<nDim;++i)
    {
      double x,y;
      a1_S3->GetPoint(i,x,y);
      //source[i] = f_s->Eval(x)/max_f;//f_s->Integral(x-0.01,x+0.01);
      int bin = h1_s3->FindBin(x);
      h1_s3->SetBinContent(bin,y);
      h1_s3->SetBinError(bin,y*0.1);
    }

  TH1F* h1_s4 = new TH1F("h1_s4","h1_s4",50,1.414,2.414);
  //double max_f = f_s->GetMaximum(1.4,2.4);
  for(int i=0;i<nDim;++i)
    {
      double x,y;
      a1_S4->GetPoint(i,x,y);
      //source[i] = f_s->Eval(x)/max_f;//f_s->Integral(x-0.01,x+0.01);
      int bin = h1_s4->FindBin(x);
      h1_s4->SetBinContent(bin,y);
      h1_s4->SetBinError(bin,y*0.1);
    }

#ifdef UNFOLD
  double Int_S3 = h1_s3->Integral("width");
  double Int_S4 = h1_s4->Integral("width");


  RooUnfoldResponse response (50, 1.414, 2.414);
  for(int i=0;i<100000;++i)
    {
      double x = gRandom->Uniform(1.450,2.340);
      double yt = f_s->Eval(x);
      double ym = a1_S3->Eval(x)/Int_S3;
      response.Fill(ym,yt);
    }
  

  RooUnfoldResponse responseS3 (50, 1.414, 2.414);
 
  for(int i=0;i<100000;++i)
    {
      double x = gRandom->Uniform(1.450,2.340);
      double yt = a1_S3->Eval(x)/Int_S3;
      double ym = a1_S4->Eval(x)/Int_S4;
      //std::cout<<"#"<<i<<" "<<x<<" "<<a1_S3->Eval(x)<<" "<<Int_S3<<" "<<y<<std::endl;
      
      responseS3.Fill(ym,yt);
    }

  TH1F* hh_responseT = (TH1F*)(responseS3.Htruth()->Clone());
  TH1F* hh_responseM = (TH1F*)(responseS3.Hmeasured()->Clone());
  TH2F* hh_responseR = (TH2F*)(responseS3.Hresponse()->Clone());

  //RooUnfoldResponse response(h1_s4,h1_s3);
  
  RooUnfoldBayes   unfold_s4 (&responseS3, h1_s4, 2);    // OR

  RooUnfoldBayes   unfold_s3 (&response, h1_s3, 2);    // OR
  //RooUnfoldSvd     unfold (&response, h1_s, 20);   // OR
  //RooUnfoldTUnfold unfold (&response, h1_s);
  
  TH1D* hReco_s4= (TH1D*) unfold_s4.Hreco();
  TH1D* hReco_s3= (TH1D*) unfold_s3.Hreco();


#endif
  // TSpectrum* spec = new TSpectrum(1);
  // char* res = spec->Deconvolution(source,response,nDim,5,2,1); 

  // std::cout<<res<<std::endl;

  // for(int i =1;i<=nDim;++i)
  //   {
  //     h1->SetBinContent(i,source[i-1]);
  //   }


  TCanvas* c1 = new TCanvas("c1","c1",500,500);
  c1->cd();
  a1_S3->Draw("a*");
  //h1_s3->Draw("same");
  //f_s->SetParameter(0,Int_S3*1.3);
  //f_s->Draw("same");

  TCanvas* c2 = new TCanvas("c2","c2",500,500);
  c2->cd();
  a1_S4->Draw("a*");
  //h1_s4->Draw("same");

#ifdef UNFOLD
  TCanvas* c3 = new TCanvas("c3","c3",500,500);
  c3->Divide(2,2);
  c3->cd(1);
  hh_responseR->Draw("colz");
  c3->cd(2);
  hh_responseT->Draw("");
  c3->cd(3);
  hh_responseM->Draw("");
  c3->cd(4);
  f_s->Draw();
  //hh_response->Draw();
  //response.Hmeasured()->Draw("");
  // a1_S3->Draw("a*");
  // h1_s3->Draw("same");

  TCanvas* c4 = new TCanvas("c4","c4",500,500);
  c4->cd();
  hReco_s3->Draw();
  hReco_s4->SetLineColor(kRed);
  hReco_s4->Draw("same");
#endif
}


void LoadResBrho(std::string namefile)
{

  std::ifstream ifs ( namefile.c_str() );
  std::string temp_line;
  std::vector<std::map<double, std::vector<double> > > energy_brho(3);
  std::vector<std::map<double,double > > brhoAll(3);

  std::getline(ifs,temp_line);

  while(std::getline(ifs,temp_line))
    {
      std::stringstream stream(temp_line);
      double energy,S,brho,sigma_brho;
     
      stream >> energy >> S >> brho >> sigma_brho;
      std::vector<double> temp_vec(2);
      temp_vec[0]=brho;
      temp_vec[1]=sigma_brho;
      //std::cout<<" S:"<<S<<" "<<energy<<" "<<eff<<" "<<eff_tot<<std::endl;
      energy_brho[S-2].insert(std::pair<double,std::vector<double> >(energy,temp_vec));
      brhoAll[S-2].insert(std::pair<double,double>(brho,sigma_brho));

    }

  TF1* f_s = new TF1("f1","[0]*TMath::Gaus(x,[1],[2],1)",.3,3.4);
  f_s->SetParameter(0,1);
  f_s->SetParameter(1,1.805088);
  f_s->SetParameter(2,0.1801583);

  TGraph* a1_S3 = new TGraph;
  a1_S3->SetNameTitle("S3_Brho","S3_Brho");
  TGraph* a1_S4 = new TGraph;
  a1_S4->SetNameTitle("S4_Brho","S4_Brho");

  TGraph* a1_S3_s = new TGraph;
  a1_S3_s->SetNameTitle("S3_BrhoSigma","S3_BrhoSigma");
  TGraph* a1_S4_s = new TGraph;
  a1_S4_s->SetNameTitle("S4_BrhoSigma","S4_BrhoSigma");

  TGraph* a1_brho = new TGraph;
  a1_brho->SetNameTitle("Brho_Sigma","Brho_Sigma");

  for(std::map<double,std::vector<double> >::iterator it_energy = energy_brho[1].begin(),it_energy_end = energy_brho[1].end();it_energy!=it_energy_end;++it_energy)
    {
      //std::map<double,double>::iterator it_eff_S2 = energy_eff[0].find(it_energy->first);
      double brho_S3 = it_energy->second[0];
      double brhoSigma_S3 = it_energy->second[1];

      //double int_f_s = f_s->Integral(it_energy->first*1.e-3-0.01,it_energy->first*1.e-3+0.01);
      if(brho_S3>1e-4)
	{
	  a1_S3->SetPoint(a1_S3->GetN(),it_energy->first*1.e-3,brho_S3);
	  a1_S3_s->SetPoint(a1_S3_s->GetN(),it_energy->first*1.e-3,brhoSigma_S3/brho_S3);
	}
    }

  for(std::map<double,std::vector<double> >::iterator it_energy = energy_brho[2].begin(),it_energy_end = energy_brho[2].end();it_energy!=it_energy_end;++it_energy)
    {
      double brho_S4 = it_energy->second[0];
      double brhoSigma_S4 = it_energy->second[1];

      //double int_f_s = f_s->Integral(it_energy->first*1.e-3-0.01,it_energy->first*1.e-3+0.01);
      if(brho_S4>1e-4)
	{
	  a1_S4->SetPoint(a1_S4->GetN(),it_energy->first*1.e-3,brho_S4);
	  a1_S4_s->SetPoint(a1_S4_s->GetN(),it_energy->first*1.e-3,brhoSigma_S4/brho_S4);
	}
    }
  for(std::map<double,double>::iterator it_brho = brhoAll[2].begin(),it_brho_end = brhoAll[2].end();it_brho!=it_brho_end;++it_brho)
    {
      double brho = it_brho->first;
      double brhoSigma = it_brho->second;

      //double int_f_s = f_s->Integral(it_energy->first*1.e-3-0.01,it_energy->first*1.e-3+0.01);
      if(brho>1e-4)
	a1_brho->SetPoint(a1_brho->GetN(),brho,brhoSigma/brho);
    }

  //TF1* f_s = new TF1("f1","TMath::Gaus(x,[0],[1],[2])",.3,3.4);
  

  TCanvas* c1 = new TCanvas("c1","c1",1000,500);
  c1->Divide(2);
  c1->cd(1);
  a1_S3->Draw("a*");
  c1->cd(2);
  a1_S3_s->Draw("a*");
  //h1_s3->Draw("same");
  //f_s->SetParameter(0,Int_S3*1.3);
  //f_s->Draw("same");

  TCanvas* c2 = new TCanvas("c2","c2",500,500);
  c2->Divide(2);
  c2->cd(1);
  a1_S4->Draw("a*");
  c2->cd(2);
  a1_S4_s->Draw("a*");
  //h1_s4->Draw("same");

  TCanvas* c3 = new TCanvas("c3","c3",500,500);
  c3->cd(1);
  a1_brho->Draw("a*");


}
