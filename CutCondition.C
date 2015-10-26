#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TCanvas.h"
#include "TMath.h"

double gauss2D(double *x, double *par) 
{
  double z1 = double((x[0]-par[1])/par[2]);
  double z2 = double((x[1]-par[3])/par[4]);
  return par[0]*exp(-0.5*(z1*z1+z2*z2));
}

int CutCondition(void)
{

  TF1* f1 = new TF1("f1","gaus(0)",-30,30);
  TF1* f2 = new TF1("f2","gaus(0)",-30,30);
  TF1* f3 = new TF1("f3","gaus(0)",-30,30);

  f1->SetParameter(0,0.75);
  f1->SetParameter(1,4.535090e-01);
  f1->SetParameter(2,4.715422e+0);
  f2->SetParameter(0,0.27);
  f2->SetParameter(1,5.122452e+00);
  f2->SetParameter(2,3.382923e0);
  f3->SetParameter(0,0.445);
  f3->SetParameter(1,1.317771e+00);
  f3->SetParameter(2,4.911264e+00);

  double tot1 = f1->Integral(-30.,30.);
  double tot2 = f2->Integral(-30.,30.);
  double tot3 = f3->Integral(-30.,30.);

  //std::vector<double> integral1;
  //std::vector<double> integral2;
  //std::vector<double> ratio;

  int size = 600;

  TGraph* g1 = new TGraph(size); 
  TGraph* g2 = new TGraph(size); 
  TGraph* g3 = new TGraph(size); 
  TGraph* g4 = new TGraph(size); 
  TGraph* g5 = new TGraph(size); 

  for(int i=0;i<size;++i)
    {
      double temp_cut = -30+i*0.1;
      double temp_int1 = f1->Integral(-30.,temp_cut);
      //integral1.push_back(temp_int1);
      double temp_int2 = f2->Integral(-30.,temp_cut);
      double temp_int3 = f3->Integral(-30.,temp_cut);
      //Integral2.push_back(temp_int2);
      //ratio.push_back(temp_int1/(temp_int1+temp_int2));
      g1->SetPoint(i,temp_cut,temp_int1/tot1);
      g2->SetPoint(i,temp_cut,1-temp_int2/tot2);
      g3->SetPoint(i,temp_cut,1-temp_int3/tot3);
      if(TMath::Abs(temp_int1+temp_int2)>0)
	{
	  g4->SetPoint(i,temp_cut,temp_int1/(temp_int1+temp_int2));
	  g5->SetPoint(i,temp_cut,temp_int1*temp_int2/((temp_int1+temp_int2)*(temp_int1+temp_int2)));
	}
      else
	{
	  g4->SetPoint(i,temp_cut,1);
	}
    }

  TMultiGraph* mg = new TMultiGraph();
  g4->SetLineColor(3);
  g5->SetLineColor(4);
  mg->Add(g1);
  mg->Add(g2);
  mg->Add(g3);
  mg->Add(g4);
  mg->Add(g5);
  
  f2->SetLineColor(3);
  f3->SetLineColor(4);
  TCanvas* c0 = new TCanvas("c0","c0",500,500);
  c0->Divide(2,2);
  c0->cd(1);
  f1->Draw();
  c0->cd(2);
  f2->Draw();
  c0->cd(3);
  f3->Draw();
  c0->cd(4);
  f1->Draw();
  f2->Draw("same");
  f3->Draw("same");
  TCanvas* c1 = new TCanvas("c1","c1",500,500);
  c1->cd();
  g1->Draw("al");
  TCanvas* c2 = new TCanvas("c2","c2",500,500);
  c2->cd();
  g2->Draw("al");
  TCanvas* c3 = new TCanvas("c3","c3",500,500);
  c3->cd();
  g3->Draw("al");
  TCanvas* c4 = new TCanvas("c4","c4",500,500);
  c4->cd();
  g4->Draw("al");
  TCanvas* c5 = new TCanvas("c5","c5",500,500);
  c5->cd();
  g5->Draw("al");
  TCanvas* c6 = new TCanvas("c6","c6",500,500);
  c6->cd();
  mg->Draw("al");
  
  return 0;
}


int CutCondition2D(void)
{

  double iniParams1[5] = {0.75,4.535090e-01,4.715422e+0,-2.922457e-03,5.057761e-01};
  TF2 * f1 = new TF2("f1",gauss2D,-30,30,-10,10, 5);
  f1->SetParameters(iniParams1);

  double iniParams2[5] = {  0.27,5.122452e+00,3.382923e0,3.014174e-03,4.708431e-01};
  TF2 * f2 = new TF2("f2",gauss2D,-30,30,-10,10, 5);
  f2->SetParameters(iniParams2);

  TCanvas* c1 = new TCanvas();
  c1->cd();
  f1->Draw("cont2");
  f2->Draw("cont2 same");









}
