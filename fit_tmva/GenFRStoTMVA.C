#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TRandom3.h"
#include <cassert>

#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

struct DataKinematics {

  float x;
  float a;
  float y;
  float b;
  float E;
  float time;
  float mass;
  float Z;
  float tof;
  float Brho;
};

struct DataIon{
  int ion;
  DataKinematics initialS;
  DataKinematics finalS;
};

int GenTMVATree(TString name_in,TString name_out)//, Long64_t N_Split, Long64_t Nb_ions, int lastN)
{


  TFile* f_in = new TFile(name_in);
  TTree* t_in = (TTree*) f_in->Get("h1");
  if(t_in==0)
    {
      cout<<"TTree empty !"<<name_in<<endl;
      f_in->ls();
      return -1;
    }

  TTreeReader reader(t_in);

  TTreeReaderValue<Int_t> N_In(reader,"N");
  TTreeReaderArray<Float_t> X_In(reader, "X");
  TTreeReaderArray<Float_t> A_In(reader, "A");
  TTreeReaderArray<Float_t> Y_In(reader, "Y");
  TTreeReaderArray<Float_t> B_In(reader, "B");
  TTreeReaderArray<Float_t> Energy_In(reader, "Energy");
  TTreeReaderArray<Float_t> Time_In(reader, "Time");
  TTreeReaderArray<Float_t> Mass_In(reader, "Mass");
  TTreeReaderArray<Float_t> Z_In(reader, "Z");
  TTreeReaderArray<Float_t> Tof_In(reader, "Tof");
  TTreeReaderArray<Float_t> Brho_In(reader, "Brho");

  const Long64_t Entries_In = reader.GetEntries(true);

  std::cout<<" Entries Input:"<<Entries_In<<std::endl;

  reader.Next();
  float mass_test = Mass_In[0];
  float z_test = Z_In[0];
  
  //std::vector<Long64_t> N_Splits;
  
  Long64_t nSplit_test = 1;
  int N_FRSstage = 0;
  while(reader.Next())
    {
      if(N_FRSstage<*N_In)
	N_FRSstage = *N_In;
      if(TMath::Abs(mass_test-Mass_In[0])<1e-1 && TMath::Abs(z_test-Z_In[0])<1e-1)
	++nSplit_test;
      // else
      // 	{
      // 	  N_Splits.emplace_back(nSplit_test);
      // 	  nSplit_test = 1;
      // 	  ++Nion_test;
      // 	  mass_test = Mass_In[0];
      // 	  z_test = Z_In[0];
      // 	}
    }
  
  int Nion_test = Entries_In / nSplit_test;
  
   //assert(N_Splits.size() == Nion_test);
  
  Long64_t N_Split = nSplit_test;
  // for(auto NS : N_Splits)
  //   if(NS != N_Split)
  //     {
  // 	std::cout<<"E> N_Split are different !"<<NS<<" "<<N_Split<<"\n";
  // 	return -1;
  //     }
  Long64_t Nb_ions = Nion_test; 
  int lastN = N_FRSstage-1;
  
  assert(Entries_In == N_Split*Nb_ions);//, "Total Entries is not multiple of N_Split and Nb_ions");
  
  //std::vector< std::vector< DataIon> > Inputs ( N_Split, std::vector<DataIon>(Nb_ions));
  
  
  TFile  *fFile = new TFile(name_out,"RECREATE");  
  TTree  *fDataT;// = new TTree("MonteCarlo", "Filtered Monte Carlo Events");
  Float_t x_i, a_i, y_i, b_i, Energy_i, Time_i, Tof_i,Brho_i,
          x_f, a_f, y_f, b_f, Energy_f, Time_f, Tof_f, Brho_f, Mass, Z;
     
  Long64_t Nevent;
  Int_t Nion;
  //Int_t IsPion, type;
  
  fDataT = new TTree("Data_Tree", "Filtered data Events");

  fDataT->Branch("x_i", &x_i, "x_i/F");
  fDataT->Branch("a_i", &a_i, "a_i/F");
  fDataT->Branch("y_i", &y_i, "y_i/F");
  fDataT->Branch("b_i", &b_i, "b_i/F");
  fDataT->Branch("Energy_i",  &Energy_i,  "Energy_i/F");
  fDataT->Branch("Time_i",  &Time_i,  "Time_i/F");
  fDataT->Branch("Tof_i",  &Tof_i,  "Tof_i/F");
  fDataT->Branch("Brho_i",  &Brho_i,  "Brho_i/F");
  fDataT->Branch("x_f", &x_f, "x_f/F");
  fDataT->Branch("a_f", &a_f, "a_f/F");
  fDataT->Branch("y_f", &y_f, "y_f/F");
  fDataT->Branch("b_f", &b_f, "b_f/F");
  fDataT->Branch("Energy_f",  &Energy_f,  "Energy_f/F");
  fDataT->Branch("Time_f",  &Time_f,  "Time_f/F");
  fDataT->Branch("Tof_f",  &Tof_f,  "Tof_f/F");
  fDataT->Branch("Brho_f",  &Brho_f,  "Brho_f/F");

  fDataT->Branch("Mass", &Mass, "Mass/F");
  fDataT->Branch("Z", &Z, "Z/F");
  fDataT->Branch("Nevent", &Nevent, "NEvent/I");
  fDataT->Branch("Nion", &Nion, "Nion/I");

  Long64_t total_nentries= N_Split;

  cout<<"N event :"<<total_nentries<<endl;
  
  Long64_t timing=0; 


  for(Long64_t i = 0;i<total_nentries;++i)
    {
      //t_in->GetEntry(i);
      if(i%100000==0)
	std::cout<<"Processing Event#"<<i<<" | "<<(double)i/(double)(total_nentries)*100<<" %"<<std::endl;
      
      if((int)((double)i/(double)(total_nentries)*10)==timing)
	{
	  std::cout<<"Progress :"<<(int)((double)i/(double)(total_nentries)*100.)<<" %"<<std::endl;
	  ++timing; 
	}

      //std::vector<DataIon> Ions;

      for(Long64_t k = 0; k<Nb_ions;++k)
	{
	  auto status = reader.SetEntry(i+k*N_Split);
	  if(status == 0)
	    {
	      
	      if((*N_In)-1!=lastN)
		continue;
	      else if((*N_In)-1 == lastN)
		{
		  //Ions.emplace_back({k,
		  //	{N.[0],X.[0],A.[0],Y.[0],B.[0],Energy.[0],Time.[0],Mass.[0],Z.[0],Tof.[0],Brho.[0]},
		  //	  {N.[lastN],X.[lastN],A.[lastN],Y.[lastN],B.[lastN],Energy.[lastN],Time.[lastN],Mass.[lastN],Z.[lastN],Tof.[lastN],Brho.[lastN]}	});

		  Nion = k;
		  x_i = X_In[0];
		  a_i = A_In[0];
		  y_i = Y_In[0];
		  b_i = B_In[0];
		  Energy_i = Energy_In[0];
		  Time_i = Time_In[0];
		  Tof_i = Tof_In[0];
		  Brho_i =  Brho_In[0];

		  x_f = X_In[lastN];
		  a_f = A_In[lastN];
		  y_f = Y_In[lastN];
		  b_f = B_In[lastN];
		  Energy_f = Energy_In[lastN];
		  Time_f = Time_In[lastN];
		  Tof_f = Tof_In[lastN];
		  Brho_f =  Brho_In[lastN];


		  Mass = Mass_In[0];
		  Z = Z_In[0];
		  Nevent = i;
		  fDataT->Fill();
      
		}
	    }	    
	}
    }

  fFile->Write();
  
  return 0;

}
