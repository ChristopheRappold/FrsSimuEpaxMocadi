/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVARegressionApplication                                          *
 *                                                                                *
 * This macro provides a simple example on how to use the trained regression MVAs *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <list>
#include <iostream>
#include <map>
#include <string>
#include <algorithm>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TH1F.h"
#include "TH2F.h"


#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

using namespace TMVA;

void TMVARegressionApplication(TString nameIn, float Mass_sel, float Z_sel, Long64_t Nevent = 0, TString myMethodList = "" ) 
{
   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   Use["PDERS"]           = 0;
   //Use["PDEFoam"]         = 1; 
   Use["KNN"]             = 1;
   // 
   // --- Linear Discriminant Analysis
   Use["LD"]		        = 1;
   // 
   // // --- Function Discriminant analysis
   // Use["FDA_GA"]          = 1;
   // Use["FDA_MC"]          = 0;
   // Use["FDA_MT"]          = 0;
   // Use["FDA_GAMT"]        = 0;
   // 
   // --- Neural Network
   Use["MLP"]             = 1; 
   // 
   // --- Support Vector Machine 
   Use["SVM"]             = 0;
   // 
   // --- Boosted Decision Trees
   Use["BDT"]             = 0;
   Use["BDTG"]            = 1;
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVARegressionApplication" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (auto it = Use.begin(); it != Use.end(); ++it)
	it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (auto it = Use.begin(); it != Use.end(); ++it)
	      std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   for(auto it = Use.begin(), ite = Use.end(); it != ite;)
     {
       if(it->second == 0)
	 it = Use.erase(it);
       else
	 ++it;
     }


   // --------------------------------------------------------------------------------------------------

   // --- Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   Float_t x_i, a_i, y_i, b_i, Energy_i, Time_i, Tof_i,Brho_i, x_f, a_f, y_f, b_f, Energy_f, Time_f, Tof_f, Brho_f, Mass, Z;
   Int_t Nion;

   reader->AddVariable( "x_i",&x_i); 
   reader->AddVariable( "a_i",&a_i); 
   reader->AddVariable( "y_i",&y_i); 
   reader->AddVariable( "b_i",&b_i); 
   reader->AddVariable( "x_f",&x_f); 
   reader->AddVariable( "a_f",&a_f); 
   reader->AddVariable( "y_f",&y_f); 
   reader->AddVariable( "b_f",&b_f); 

   // Spectator variables declared in the training have to be added to the reader, too
   reader->AddSpectator( "Mass", &Mass);
   reader->AddSpectator( "Z", &Z);
   reader->AddSpectator( "Nion",&Nion);

   // --- Book the MVA methods

   TString dir    = "WeightsReg/";
   TString prefix = "TMVARegression";

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); ++it) {
      if (it->second) {
         TString methodName = it->first + " method";
         TString weightfile = dir + prefix + "_" + TString(it->first) + ".weights.xml";
         reader->BookMVA( methodName, weightfile ); 
      }
   }
   
   // Book output histograms
   std::map<std::string, TH1F*> hists ;
   std::map<std::string, TH2F*> hists2D ;
   std::map<std::string, TH2F*> hists2D_Diff ;
   std::map<std::string, TH2F*> hists2D_Pull ;
   std::map<std::string, TH2F*> hists2D_Power ;
   std::string methods (" method");
   std::list<std::string> UsedMethod;
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); ++it)
     {
       if (it->second)
	 {
	   TString name_temp = it->first;
	   std::string name_method = it->first;
	   name_method += methods;
	   name_temp += "_Value";
	   TH1F* h = new TH1F( name_temp, name_temp, 100, 13, 15 );
	   hists.insert(std::make_pair(name_method, h));       

	   TString name_temp1 = it->first;
	   name_temp1 += "_FitVsReal";
	   TH2F* h2_0 = new TH2F(name_temp1,name_temp1,100,13,15,100,13,15);
	   hists2D.insert(std::make_pair(name_method, h2_0));       
	   
	   TString name_temp2 = it->first;
	   name_temp2 += "_Diff";
	   TH2F* h2_1 = new TH2F(name_temp2,name_temp2,100,13,15,1000,-0.1,0.1);
	   hists2D_Diff.insert(std::make_pair(name_method, h2_1));       

	   TString name_temp3 = it->first;
	   name_temp3 += "_Pull";
	   TH2F* h3 = new TH2F(name_temp3,name_temp3,100,13,15,1000,-0.1,0.1);
	   hists2D_Pull.insert(std::make_pair(name_method, h3));       

	   TString name_temp4 = it->first;
	   name_temp4 += "_Power";
	   TH2F* h4 = new TH2F(name_temp4,name_temp4,100,13,15,1000,0,4000);
	   hists2D_Power.insert(std::make_pair(name_method, h4));       

	   UsedMethod.emplace_back(name_method);
	 }
     }
   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //   
   TFile *input = TFile::Open(nameIn);

   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVARegressionApp        : Using input file: " << input->GetName() << std::endl;

   // --- Event loop

   // Prepare the tree
   // - here the variable names have to corresponds to your tree
   // - you can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   TTree* theTree = (TTree*)input->Get("Data_Tree");
   std::cout << "--- Select signal sample" << std::endl;

   theTree->SetBranchAddress( "x_i",&x_i); 
   theTree->SetBranchAddress( "a_i",&a_i); 
   theTree->SetBranchAddress( "y_i",&y_i); 
   theTree->SetBranchAddress( "b_i",&b_i); 
   theTree->SetBranchAddress( "x_f",&x_f); 
   theTree->SetBranchAddress( "a_f",&a_f); 
   theTree->SetBranchAddress( "y_f",&y_f); 
   theTree->SetBranchAddress( "b_f",&b_f); 
   theTree->SetBranchAddress( "Brho_f",&Brho_f);

   theTree->SetBranchAddress( "Mass",&Mass);
   theTree->SetBranchAddress( "Z",&Z);
   
   Long64_t TotEntries = Nevent == 0 ? theTree->GetEntries() : Nevent;
   
   std::cout << "--- Processing: " << TotEntries<<" / "<< theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   for (Long64_t ievt=0; ievt<TotEntries;ievt++)
     {
       
       if (ievt%1000 == 0)
	 {
	   std::cout << "--- ... Processing event: " << ievt << std::endl;
	 }

       theTree->GetEntry(ievt);
       if(TMath::Abs(Mass-Mass_sel)>1e-1 && TMath::Abs(Z-Z_sel)>1e-1)
	 continue;
       // Retrieve the MVA target values (regression outputs) and fill into histograms
       // NOTE: EvaluateRegression(..) returns a vector for multi-target regression

       for(const auto& it_use : UsedMethod)
	 {

	   TString title = it_use;
	   Float_t val = (reader->EvaluateRegression( title ))[0];

	   hists[it_use]->Fill(val);
	   hists2D[it_use]->Fill(Brho_f,val);
	   hists2D_Diff[it_use]->Fill(Brho_f,(Brho_f-val));
	   hists2D_Pull[it_use]->Fill(Brho_f,TMath::Abs(Brho_f-val)/Brho_f);
	   if(TMath::Abs(Brho_f-val)>1e-6)
	     hists2D_Power[it_use]->Fill(Brho_f,Brho_f/TMath::Abs(Brho_f-val));
	 }
     }
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   // --- Write histograms

   TFile *target  = new TFile( "TMVARegApp2.root","RECREATE" );
   target->cd();
   for(auto& it : hists)
     it.second->Write();
   for(auto& it : hists2D)
     it.second->Write();
   for(auto& it : hists2D_Diff)
     it.second->Write();
   for(auto& it : hists2D_Pull)
     it.second->Write();
   for(auto& it : hists2D_Power)
     it.second->Write();
   
   target->Close();

   std::cout << "--- Created root file: \"" << target->GetName() 
             << "\" containing the MVA output histograms" << std::endl;
  
   delete reader;
    
   std::cout << "==> TMVARegressionApplication is done!" << std::endl << std::endl;
}

// int main( int argc, char** argv )
// {
//    // Select methods (don't look at this code - not of interest)
//    TString methodList; 
//    for (int i=1; i<argc; i++) {
//       TString regMethod(argv[i]);
//       if(regMethod=="-b" || regMethod=="--batch") continue;
//       if (!methodList.IsNull()) methodList += TString(","); 
//       methodList += regMethod;
//    }
//    TMVARegressionApplication(methodList);
//    return 0;
// }
