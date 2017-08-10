#include "TH1F.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraphAsymmErrors.h"

#include <fstream>
#include <iostream>
#include <vector>

void AnaSysEnergy(std::string nameList,int min, int max)
{
  std::ifstream List(nameList.c_str());
  std::string infiles;
  int nb_file = 0;
  
  TGraphAsymmErrors* graph_sum = new TGraphAsymmErrors();
  graph_sum->SetNameTitle("Summary_maxEff","Summary_maxEff");
  std::vector<TGraphAsymmErrors*> vec_graph;
  int index = 0;
  while(std::getline(List, infiles))
    {
      std::cout << infiles << std::endl;
      TFile* f_in = new TFile(infiles.c_str());
      f_in->cd();
      TTree* tree = (TTree*)f_in->Get("T");
      TString name1("brho>>h1(60,");
      name1+=min;
      name1+=",";
      name1+=max;
      name1+=")";
      TString name2("brho>>h2(60,");
      name2+=min;
      name2+=",";
      name2+=max;
      name2+=")";
      tree->Draw(name1.Data(),"n==1","goff");
      tree->Draw(name2.Data(),"n==2","goff");

      TH1F* h1 = (TH1F*)gDirectory->Get("h1");
      TString nameHistAll("HistAll_");
      nameHistAll+=index;
      h1->SetNameTitle(nameHistAll,nameHistAll);

      TH1F* h2 = (TH1F*)gDirectory->Get("h2");
      TString nameHistAcc("HistAcc_");
      nameHistAcc+=index;
      h2->SetNameTitle(nameHistAcc,nameHistAcc);

      
      TH1F* h11 = (TH1F*)h1->Clone();
      TH1F* h22 = (TH1F*)h2->Clone();
      h11->Sumw2();
      h22->Sumw2();
      
      TGraphAsymmErrors* g_acc = new TGraphAsymmErrors(h22,h11);
      TString nameAcc("Eff_");
      nameAcc+=index;
      g_acc->SetNameTitle(nameAcc,nameAcc);
      vec_graph.push_back(g_acc);

      double x_max=0,y_max=-1, x_errL_max=0, x_errH_max=0, y_errL_max=0, y_errH_max=0;
      for(int bin_i = 0; bin_i<g_acc->GetN();++bin_i)
	{
	  double x=0.,y=0.;
	  g_acc->GetPoint(bin_i,x,y);
	  if(y_max<y)
	    {
	      x_max = x;
	      y_max = y;
	      x_errL_max = g_acc->GetErrorXlow(bin_i);
	      x_errH_max = g_acc->GetErrorXhigh(bin_i);
	      y_errL_max = g_acc->GetErrorYlow(bin_i);
	      y_errH_max = g_acc->GetErrorYhigh(bin_i);
	    }
	}
      if(y_max>0.)
	{
	  graph_sum->SetPoint(index,x_max,y_max);
	  graph_sum->SetPointError(index,x_errL_max,x_errH_max,y_errL_max,y_errH_max);
	}
      ++index;
      f_in->Close();
    }
  
  std::cout << " Loaded " << nb_file << " files " << std::endl;
  TFile* f_out = new TFile("AllResEff.root","RECREATE");
  f_out->cd();
  graph_sum->Write();
  for(size_t i = 0; i<vec_graph.size();++i)
    {
      vec_graph[i]->Write();
    }
  f_out->Close();
}

