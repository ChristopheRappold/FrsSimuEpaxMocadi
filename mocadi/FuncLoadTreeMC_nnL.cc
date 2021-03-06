#include "TFile.h"
#include "TMath.h"
#include "TTree.h"
#include "Riostream.h"
#include "TVector3.h"
#include "TLorentzVector.h"
//#include "TROOT.h"
//#include "TSystem.h"
// #include "TH1.h"
// #include "TH2.h"

#include "MCAnaEventG4Sol.hh"
#include "TMcParticle.hh"
#include "TMcHit.hh"
#include "THyphiTrack_v4.hh"

static TFile* f = 0;
static TTree* t = 0;
static MCAnaEventG4Sol* event = 0;
static int id_event = 0;

//int FuncLoadTree(double *in, double *out, double *dpar, char *option);

extern "C" int FuncLoadTree(double *in, double *out, double *dpar, char *option)
{
 /* in[0]=X [cm]   in[4]=energy[AMeV] in[8]=electron       in[12]=deltaE[MeV]
     in[1]=X'[mrad] in[5]=time  [us]   in[9]=nf/nsf         in[13]=reserved
     in[2]=Y [cm]   in[6]=mass  [amu]  in[10]=range[mg/c2]  in[13]=reserved
     in[3]=Y'[mrad] in[7]=z            in[11]=tof  [us]
     the element range is valid after the "stop" keyword
     the element deltaE is valid behind energy loss materials
                                          (matter, wedge etc.)
  */

  if(f==0)
    {
      f= new TFile(option);
      event = 0;
      t= dynamic_cast<TTree*>(f->Get("T"));
      t->SetCacheSize(30000000);
      t->AddBranchToCache("*");
      t->SetBranchAddress("MCAnaEventG4Sol",&event);
    }
  
  int start_ev = (int)(dpar[0]);
  if(id_event%5000==1)
    std::cout<<"current event#"<<id_event+start_ev<<"\n";

  
  t->GetEntry(id_event+start_ev);
  ++id_event;
  int d_id = -1;
  if(event->Nmc!=0)
    {
      for(int k = 0 ; k<event->Nmc ;++k)
	{
	  TMcParticle* par = dynamic_cast<TMcParticle*>(event->fMC_Particle->At(k));
	  if(par->Charge==1)
	    d_id = par->Mc_id;
	}
    }
  if(d_id==-1)
    return -1;

  int ok = -1;
  for(int k = 0 ; k<event->FMF2->GetEntries() ;++k)
    {
      TMcHit* hitStop = dynamic_cast<TMcHit*>(event->FMF2->At(k));
      
      if(hitStop->MC_id==d_id && hitStop->LayerID==47)
	{
	  ok = 1;
	  double mass_u = hitStop->MCparticle.M()/0.931494028; // amu
	  out[0] = hitStop->MCHit.X(); // cm
	  out[1] = TMath::ATan(hitStop->MCparticle.X()/hitStop->MCparticle.Z())*1000.; // mrad
	  out[2] = hitStop->MCHit.Y(); // cm
	  out[3] = TMath::ATan(hitStop->MCparticle.Y()/hitStop->MCparticle.Z())*1000.; // mrad
	  out[4] = (hitStop->MCparticle.E()-hitStop->MCparticle.M())*1000./mass_u; // AMeV
	  out[5] = 0. ; // micros
	  out[6] = mass_u;
	  out[7] = 1; // charge
	  out[8] = 0; // electron
	  out[9] = in[9]; // nf/nsf
	  out[10] = in[10];
	  out[11] = 0.;
	  out[12] = in[12];
	  out[13] = in[13];
	  out[14] = in[14];
	  //point->SetPoint(IdN,hitStop->MCHit.X(),hitStop->MCHit.Y(),hitStop->MCHit.Z());
	  //++IdN;
	  //std::cout<<" STOPx:"<<hitStop->MCHit.X()<<", STOPy"<<hitStop->MCHit.Y()<<" STOPz:"<<hitStop->MCHit.Z()<<" Charge:"<<hitStop->Charge<<std::endl;
	  //}
	}
    }
  if(ok==1)
    return 0;
  else
    return -1;



}
