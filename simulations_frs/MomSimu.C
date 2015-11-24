// #include "../../FullReco_Mainz/src/EventG4/THypHi_Par.h"
// #include "../../FullReco_Mainz/src/EventG4/THypHi_FiberPlain.h"
// #include "../../FullReco_Mainz/src/EventG4/THypHi_Event.h"
// #include "../../FullReco_Mainz/src/EventG4/TUTracker_Event.h"
#include "../../FullReco_Mainz/src/Ana_Event/MCAnaEvent.hh"
#include "../../FullReco_Mainz/src/Ana_Event/TMcParticle.hh"
#include "TMath.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TCanvas.h"



void MomSimu(int index2, TString name_in, Long64_t nb=0)
{
  TF1* f0 = new TF1("f0","1./TMath::Sqrt(2.*TMath::Pi())/[1]*TMath::Exp(-0.5*(x-[0])*(x-[0])/[1]/[1])",-10,10);
  TF1* f1 = new TF1("f1","x/TMath::Sqrt(2.*TMath::Pi())/[1]*TMath::Exp(-0.5*(x-[0])*(x-[0])/[1]/[1])",-10,10);

  int index  = 3; // t ; He3 ; He4 ; d

  double min[4] = {6.98,9.1,12,4.48};
  double max[4] = {9.25,12.2,17.55,6.16};

  
  double mean_p[4] = {8.177,10.625,14.82,5.315};
  double sigma_p[4] = {0.286,0.4046,0.7584,0.2849};

  
  f0->SetParameter(0,mean_p[index]);
  f0->SetParameter(1,sigma_p[index]);
  f1->SetParameter(0,mean_p[index]);
  f1->SetParameter(1,sigma_p[index]);

  double temp_diff = 1.;
  double interval_max = max[index];
  double previous_max = max[index];
  int i = 0;
  
  while(TMath::Abs(temp_diff-0.01)>1e-5)
    {
      double temp_mean = previous_max-(previous_max-min[index])*(TMath::Abs(temp_diff-0.01));
      double temp_max = temp_mean;
      //cout<<"---"<<endl;
      temp_diff = TMath::Sqrt(f0->Variance(min[index],temp_max))/f0->Mean(min[index],temp_max);
      if(temp_diff>0.01)
	{
	  previous_max = interval_max;
	  interval_max = temp_max;
	  //cout<<" diff big "<<temp_diff<<endl;
	}
      if(temp_diff<0.01)
	{
	  //cout<<" diff small"<<temp_diff<<" "<<temp_max<<" | ";
	  temp_max = 0.5*(temp_max+previous_max);
	  temp_diff = TMath::Sqrt(f0->Variance(min[index],temp_max))/f0->Mean(min[index],temp_max);
	  if(temp_diff>0.01)
	    previous_max = temp_max;
	  cout<<temp_max<<" "<<temp_diff<<endl;
	}

      //cout<<"i ["<<i<<"] "<<min<<" "<<temp_max<<" "<<previous_max<<" "<<temp_diff<<endl;
      ++i;
      interval_max = temp_max;
      if(i>1000)
	break;
    }
  
  cout<<" last :"<<interval_max<<" "<<TMath::Sqrt(f0->Variance(min[index],interval_max))/f0->Mean(min[index],interval_max)<<endl;
  


  TFile* f= new TFile(name_in);

  //THypHi_Event* event = 0;
  //TTree* t= (TTree*)f->Get("HypHiMC_output_TREE");
  MCAnaEvent* event = 0;
  TTree* t= (TTree*)f->Get("T");
  t->SetCacheSize(30000000);
  t->AddBranchToCache("*");

  //t->SetBranchAddress("HypHi_Event",&event);
  t->SetBranchAddress("MCAnaEvent",&event);
  

  double last_max = 0;
  double last_min = min[index];//mean_p[index]-2*sigma_p[index];//min[index]
  //double last_min = mean_p[index]-2*sigma_p[index];//min[index]
  i =0;
  std::vector< std::vector<double> > list_bondary;
  double bondary_max = max[index];//mean_p[index]+2*sigma_p[index];//max[index];
  //double bondary_max = mean_p[index]+2*sigma_p[index];//max[index];
  
  while(last_max<bondary_max)
    {
      std::vector<double> temp_vec(2,0.);
      //last_max = 1.005/0.995 * last_min;
      last_max = 1.01/0.99 * last_min;
      temp_vec[0]=last_min;
      temp_vec[1]=last_max;
      list_bondary.push_back(temp_vec);
      cout<<" ["<<last_min<<" "<<last_max<<"]"<<endl;
      last_min = last_max - 0.25*(last_max-last_min);
      ++i;

      if(i>100)
	break;
    }

  cout<<" nb interval :"<<i<<endl;

  TH1F* h1[list_bondary.size()];
  for(unsigned int k=0;k<list_bondary.size();++k)
    {
      TString name_temp("Ranged_");
      name_temp+=k;
      h1[k] = new TH1F(name_temp,name_temp,200,0.,6.7);
    }
  Long64_t Nevents = nb != 0 ? nb : t->GetEntries();
  cout<<"N :"<<Nevents<<endl;
  
  
  for(Long64_t ie = 0; ie<Nevents;++ie)
    {
      if(ie%10000==1)
	cout<<"current event#"<<ie<<endl;
      t->GetEntry(ie);
      
      //for(unsigned int k = 0;k<event->BeamNames.size();++k)
      if(event->Nmc!=0)
	{
	  for(int k = 0 ; k<event->Nmc ;++k)
	    {
	      TMcParticle* par = dynamic_cast<TMcParticle*>(event->fMC_Particle->At(k));
	      //if(event->BeamNames[k]=="triton")
	      if(par->Pdg == 10000) 
		{
		  //double mom = event->BeamMomentums_Z[k] *1e-3;
		  double mom  = par->MomMass.Pz();
		  for(unsigned int j = 0; j<list_bondary.size();++j)
		    {
		      if(list_bondary[j][0]<mom && mom<list_bondary[j][1])
			{
			  h1[j]->Fill(mom);
			}
		      
		    }
		  
		}
	    }
	}
    }

  TCanvas* c1 = new TCanvas("c1","c1",500,500);
  c1->cd();
  h1[0]->Draw();
  int color[4] = {1,2,3,4};
  for(unsigned int k=1;k<list_bondary.size();++k)
    {
      h1[k]->SetLineColor(color[k%4]);
      h1[k]->Draw("same");
    }


}
