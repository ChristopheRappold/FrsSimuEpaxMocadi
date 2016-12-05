#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TRandom3.h"

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

int GEnTMVATree(TString name_in,TString name_out, Long64_t N_Split, Long64_t Nb_ions, int lastN)
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

  const auto Entries_In = reader.GetEntries(true);

  std::cout<<" Entries Input:"<<Entries_In<<std::endl;
  
  static_assert(Entries_In = N_Split*Nb_ions);
  
  //std::vector< std::vector< DataIon> > Inputs ( N_Split, std::vector<DataIon>(Nb_ions));

  while(reader.Next()) 
    {
      auto iEntry = reader.GetCurrentEntry();
      if(static_cast<int>(static_cast<double>(iEntry)/static_cast<double>(Entries)*10)==timing)
	{
	  std::cout<<"Processing :"<<timing*10<<"% \n";
	  ++timing;
	}

      // int Nion = iEntry / N_Split;
      // int Nevent = iEntry % N_Split;

      // if(N-1 == lastN)
      // 	{
      // 	  Inputs[Nevent][Nion] = {Nion,
      // 				  {N.begin(),X.begin(),A.begin(),Y.begin(),B.begin(),Energy.begin(),Time.begin(),Mass.begin(),Z.begin(),Tof.begin(),Brho.begin()},
      // 				  {N.last(),X.last(),A.last(),Y.last(),B.last(),Energy.last(),Time.last(),Mass.last(),Z.last(),Tof.last(),Brho.last()} };

	  
      // 	}
      // else 
      // 	  Inputs[Nevent][Nion].ion =-1;


    }


  

  
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
	      
	      if(N-1!=Nlast)
		continue;
	      else if(N-1 == Nlast)
		{
		  //Ions.emplace_back({k,
		  //	{N.[0],X.[0],A.[0],Y.[0],B.[0],Energy.[0],Time.[0],Mass.[0],Z.[0],Tof.[0],Brho.[0]},
		  //	  {N.[Nlast],X.[Nlast],A.[Nlast],Y.[Nlast],B.[Nlast],Energy.[Nlast],Time.[Nlast],Mass.[Nlast],Z.[Nlast],Tof.[Nlast],Brho.[Nlast]}	});

		  Nion = k;
		  x_i = X_In[0];
		  a_i = A_In[0];
		  y_i = Y_In[0];
		  b_i = B_In[0];
		  Energy_i = Energy_In[0];
		  Time_i = Time_In[0];
		  Tof_i = Tof_In[0];
		  Brho_i =  Brho_In[0];

		  x_f = X_In[Nlast];
		  a_f = A_In[Nlast];
		  y_f = Y_In[Nlast];
		  b_f = B_In[Nlast];
		  Energy_f = Energy_In[Nlast];
		  Time_f = Time_In[Nlast];
		  Tof_f = Tof_In[Nlast];
		  Brho_f =  Brho_In[Nlast];


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
