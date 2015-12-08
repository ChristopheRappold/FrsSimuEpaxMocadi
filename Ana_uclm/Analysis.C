#include "Riostream.h"
#include "TGraph.h"
#include "TCanvas.h"

#include <iostream>
#include <fstream>
#include <sstream>
//#include <unordered_map>
//#include <unordered_set>
#include <set>
#include <map>
#include <vector>
#include <list>
#include <string>
#include <cctype>
#include <algorithm>
#include <functional>
#include <iterator>
//#include <tuple>

#include "TMath.h"
#include "TString.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TMarker.h"
//#include "TRandom3.h"

#include "Analysis.h"

static const std::string ElName2[111] = {
  "n","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg",
  "Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr",
  "Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
  "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
  "In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
  "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf",
  "Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po",
  "At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm",
  "Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs",
  "Mt","Ds" };



int AnalysisHyp(std::string namefileEpax, std::string namefileHyp, std::string Hyp, int type = 0, std::string ICtype = "I-C", std::string Mtype = "A B AB", int Ntype = -1, int NtypeComp = 1, double minHist = 0, int Nbin = 20)
{
  const std::set<std::string> StableBeam = {"H1","H2","He3","He4","Li6","Li7","Be9","B10","B11","C12","C13","N14","N15","O16","O17","O18","F19","Ne20","Ne21","Ne22","Na23","Mg24","Mg25","Mg26","Al27","Si28","Si29","Si30","P31","S32","S33","S34","S36","Cl35","Cl37","Ar36","Ar38","Ar40"};
  
  const std::map<int,int> StableBeam2 = {{1,1},{1,2},{2,3},{2,4},{3,6},{3,7},{4,9},{5,10},{5,11},{6,12},{6,13},{7,14},{7,15},{8,16},{8,17},{8,18},{9,19},{10,20},{10,21},{10,22},{11,23},{12,24},{12,25},{12,26},{13,27},{14,28},{14,29},{14,30},{15,31},{16,32},{16,33},{16,34},{16,36},{17,35},{17,37},{18,36},{18,38},{18,40}};

  std::set<std::string> RejectedBeam;
  for(auto& it_beam : StableBeam2 )
    {
      std::string nameF1("");
      nameF1+= ElName2[it_beam.first];
      nameF1+=std::to_string(it_beam.second+1);
      //std::string SnameF(nameF.Data());
      std::string nameF2("");
      nameF2+= ElName2[it_beam.first];
      nameF2+=std::to_string(it_beam.second-1);

      RejectedBeam.emplace(nameF1);
      RejectedBeam.emplace(nameF2);
    }
  std::map<std::string,double> SecondTargetDensity;
  SecondTargetDensity["C12"] = 2.21/12;
  SecondTargetDensity["Be9"] = 1.85/9;
    
  double PrimaryBeamInt = 1e10;

  std::ifstream ifs ( namefileEpax.c_str() );

  //std::string temp_line;
  int Af,Zf,Ap,Zp,At,Zt;
  double TargetDensity,xs,ProdRate;
  int Sister1,Sister2,Sister3,ZZ;
  double CmTarget, gCmTarget;
  double SurvivalRate, MeanEnergy, SigmaEnergy, Transmission;
  double prod1;
  double prod2;
  double prod3;

  double prod4;
  
  std::map<std::string,tuple_graph > BeamTarget;

  std::vector<double> intensity;
  std::vector<double> contamination;
  std::vector<double> targetLength;
  std::set<int> targetLength2;
  std::vector<double> EnergyM;
  std::vector<double> Survival;
  std::vector<double> Trans;
  std::vector<double> IntCon;

  std::set<std::string> FragAvailable;
  
  std::string temp_line;
  std::getline(ifs,temp_line);

  std::vector<double> MinMaxCheck = {1.e100,-1.};
  std::vector<std::string> NameCheck(2);
  
  while(std::getline(ifs,temp_line))
    {
      std::stringstream stream(temp_line);
      stream >> Af >> Zf >> Ap >> Zp >> At >> Zt >> TargetDensity >> xs >> ProdRate >> Sister1 >> Sister2 >> Sister3 >> ZZ
	     >> CmTarget >> SurvivalRate >> MeanEnergy >> SigmaEnergy >> Transmission >> prod1 >> prod2 >> prod3 >> prod4 ; 

      if(MeanEnergy<700.)
	continue;
      if(At == 2 && Zt == 1)
	continue;
      // if(MeanEnergy-3*SigmaEnergy < 1600.)
      // 	continue;
      if(At == 6 && Zt == 3)
	continue;
      intensity.push_back(PrimaryBeamInt*prod1*(1.-SurvivalRate));
      contamination.push_back((prod2+prod3+prod4)*(1.-SurvivalRate)*PrimaryBeamInt);
      targetLength.push_back(ProdRate*CmTarget);
      
      IntCon.push_back(PrimaryBeamInt*(prod1-prod2+prod3+prod4)*(1.-SurvivalRate));
      EnergyM.push_back(MeanEnergy);
      Survival.push_back(1.-SurvivalRate);
      Trans.push_back(Transmission);
      
      std::vector<int> idsBT ;
      idsBT.push_back(Ap);
      idsBT.push_back(Zp);
      idsBT.push_back(At);
      idsBT.push_back(Zt);
      idsBT.push_back(Af);
      idsBT.push_back(Zf);
      TString nameF("");
      nameF+= ElName2[Zf];
      nameF+=Af;
      TString nameP("");
      nameP+= ElName2[Zp];
      nameP+=Ap;
      TString nameT(""); 
      nameT+= ElName2[Zt];
      nameT+=At;

      std::string SnameF(nameF.Data());
      std::string SnameP(nameP.Data());
      std::string SnameT(nameT.Data());
      std::string SnameAll(SnameF);
      SnameAll+=SnameP;
      SnameAll+=SnameT;
      std::string SnameBT(SnameP);
      SnameBT+="+";
      SnameBT+=SnameT;

      FragAvailable.emplace(SnameF);
      
      std::map<std::string,tuple_graph>::iterator it_FBT = BeamTarget.find(SnameAll);
      if(it_FBT==BeamTarget.end())
	{
	  tuple_graph temp(SnameBT,SnameF,idsBT,CmTarget,ProdRate*CmTarget,PrimaryBeamInt*prod1*(1.-SurvivalRate),PrimaryBeamInt*prod2*(1.-SurvivalRate),PrimaryBeamInt*prod3*(1.-SurvivalRate),PrimaryBeamInt*prod4*(1.-SurvivalRate),MeanEnergy,SigmaEnergy,SurvivalRate);
	  BeamTarget.insert(std::pair<std::string,tuple_graph>(SnameAll,temp));
	}
      else
	{
	  it_FBT->second.AddValue(CmTarget,ProdRate*CmTarget,PrimaryBeamInt*prod1*(1.-SurvivalRate),PrimaryBeamInt*prod2*(1.-SurvivalRate),PrimaryBeamInt*prod3*(1.-SurvivalRate),PrimaryBeamInt*prod4*(1.-SurvivalRate),MeanEnergy,SigmaEnergy,SurvivalRate);
	}

      if(IntCon.back() < MinMaxCheck[0])
	{
	  MinMaxCheck[0] = IntCon.back();
	  NameCheck[0] = SnameAll;
	}
      if(IntCon.back() > MinMaxCheck[1])
	{
	  MinMaxCheck[1] = IntCon.back();
	  NameCheck[1] = SnameAll;
	}
      
    }



  
  std::ifstream ifsHyp (namefileHyp.c_str());
  std::set<nameBTtuple,CompBTTuple> hyp_prodAll;

  std::getline(ifsHyp,temp_line);
  double max_CX = -1 ;
  while(std::getline(ifsHyp,temp_line))
    {
      std::stringstream stream(temp_line);
      std::string hyp;
      int nb;
      stream >> hyp >> nb;
      if(hyp != Hyp)
	continue;
      std::cout<<temp_line<<std::endl;
      nb /= 2;
      for(int id = 0; id<nb ; ++id)
	{
	  std::string BT;
	  double CX1;

	  stream >> BT >> CX1;

	  std::cout<<" ---> id#"<<id<<" "<<BT<<" "<<CX1;
	  std::size_t foundPlus = BT.find("+");
	  if(foundPlus == std::string::npos )
	    {
	      std::cout<<"E> no + ! "<<id<<" "<<BT<<" "<<hyp<<std::endl;
	      return -1;
	    }
	  nameBTtuple nameBT;
	  nameBT.nameB=BT.substr(0,foundPlus);
	  nameBT.nameT=BT.substr(foundPlus+1);
	  nameBT.cross_section = CX1;
	  
	  // auto isSecondaryBeamAvailable = FragAvailable.find(nameBT.nameB);
	  // if(isSecondaryBeamAvailable != FragAvailable.end())
	  //   {
	  std::cout<<" beam:"<<nameBT.nameB<<" target:"<<nameBT.nameT<<" "<<nameBT.cross_section<<std::endl; 
	  auto ret = hyp_prodAll.emplace(nameBT);
	  if (!ret.second)
	    std::cout << "already exists in hyp_prod\n";
	  else
	    {
	      if(CX1 > max_CX )
		max_CX = CX1;
	    }
	  //   }
	  // else
	  //   std::cout<<" 2nd beam not available !"<<" beam:"<<nameBT.nameB<<" target:"<<nameBT.nameT<<" "<<nameBT.cross_section<<std::endl;
	}
    }

  auto Pred = [&] (const nameBTtuple& nameBT) -> bool
	       {
		 auto isFound = FragAvailable.find(nameBT.nameB);
		 if(isFound != FragAvailable.end())
		   return true;
		 else
		   {
		     std::cout<<" 2nd beam not available !"<<" beam:"<<nameBT.nameB<<" target:"<<nameBT.nameT<<" "<<nameBT.cross_section<<" "<<nameBT.cross_section/max_CX<<std::endl;
		     return false;
		   }
	       };
  
  std::set<nameBTtuple,CompBTTuple> hyp_prod;
  for(auto item : hyp_prodAll)
    if(Pred(item))
      {
	hyp_prod.emplace(item);
      }
  
  
  if(hyp_prod.size()==0)
    {
      std::cout<<"!> no hypernuclei selected ! "<<Hyp<<std::endl;
      return -2;
    }

  std::cout<<"Hyp required: "<<Hyp<<std::endl<<" 2nd Reaction:";
  for(auto hypBT : hyp_prod)
    {
      //if(hypBT.nameT=="C12" && hypBT.nameB=="N12")
      std::cout<<" <"<<hypBT.nameB<<"> + <"<<hypBT.nameT<<"> ("<<hypBT.cross_section<<")"<<std::endl;
    }
  std::cout<<std::endl;

  
  

  // for(std::map<std::string, tuple_graph>::const_iterator it_FBT = BeamTarget.begin(), it_FBT_end = BeamTarget.end();it_FBT!= it_FBT_end;++it_FBT)
  //   {
  //      std::cout<<"*Entry :"<<it_FBT->first<<std::endl;
  //      it_FBT->second.Print();
  //   }

  // auto it_FBTmin = BeamTarget.find(NameCheck[0]);
  // std::cout<<"* Entry min :"<<it_FBTmin->first<<std::endl;
  // it_FBTmin->second.Print();

  // auto it_FBTmax = BeamTarget.find(NameCheck[1]);
  // std::cout<<"* Entry max :"<<it_FBTmax->first<<std::endl;
  // it_FBTmax->second.Print();

  auto it_FBTtest = BeamTarget.find("N12N14Be9");
  std::cout<<"* Entry test :"<<it_FBTtest->first<<std::endl;
  it_FBTtest->second.Print();

  Compute compute;
  LambdaCrossSectionFunction StrangenessProdNorm(2000.);
  
  double Tminmax[2] = {*std::min_element(targetLength.begin(),targetLength.end()),*std::max_element(targetLength.begin(),targetLength.end())};
  double Iminmax[2] = {*std::min_element(intensity.begin(),intensity.end()),*std::max_element(intensity.begin(),intensity.end())};
  double Cminmax[2] = {*std::min_element(contamination.begin(),contamination.end()),*std::max_element(contamination.begin(),contamination.end())};
  double ICminmax[2] = {*std::min_element(IntCon.begin(),IntCon.end()),*std::max_element(IntCon.begin(),IntCon.end())};
  double Energyminmax [2] ={*std::min_element(EnergyM.begin(),EnergyM.end()),*std::max_element(EnergyM.begin(),EnergyM.end())};
  double Survivalminmax [2] ={*std::min_element(Survival.begin(),Survival.end()),*std::max_element(Survival.begin(),Survival.end())};
  double Transminmax [2] ={*std::min_element(Trans.begin(),Trans.end()),*std::max_element(Trans.begin(),Trans.end())};

  Iminmax[1] = 2.e7;
  
  double Tmeansigma [2];
  double Imeansigma [2];
  double Cmeansigma [2];
  double ICmeansigma [2];
  double Energymeansigma [2];
  double Survivalmeansigma [2];
  double Transmeansigma [2];
  
  compute.computeMeanSigma(intensity,Imeansigma);
  compute.computeMeanSigma(targetLength,Tmeansigma);
  compute.computeMeanSigma(contamination,Cmeansigma);
  compute.computeMeanSigma(IntCon,ICmeansigma);
  compute.computeMeanSigma(EnergyM, Energymeansigma);
  compute.computeMeanSigma(Survival, Survivalmeansigma);
  compute.computeMeanSigma(Trans, Transmeansigma);

  double Tquantile [3];
  double Iquantile [3];
  double Cquantile [3];
  double ICquantile [3];
  double Energyquantile [3];
  double Survivalquantile [3];
  double Transquantile [3];
  
  compute.computeQuantile(intensity,Iquantile);
  compute.computeQuantile(targetLength,Tquantile);
  compute.computeQuantile(contamination,Cquantile);
  compute.computeQuantile(IntCon,ICquantile);
  compute.computeQuantile(EnergyM, Energyquantile);
  compute.computeQuantile(Survival, Survivalquantile);
  compute.computeQuantile(Trans, Transquantile);

  
  std::vector<std::vector<double*> > RangeAll = {{Tminmax,Iminmax,Cminmax,Energyminmax,Survivalminmax,Transminmax},
						 {Tmeansigma,Imeansigma,Cmeansigma,Energymeansigma,Survivalmeansigma,Transmeansigma},
						 {Tquantile,Iquantile,Cquantile,Energyquantile,Survivalquantile,Transquantile}};

    if(ICtype == "I-C")
    {
      RangeAll[0][2] = ICminmax;
      RangeAll[1][2] = ICmeansigma;
      RangeAll[2][2] = ICquantile;
    }

  std::cout<<" mix max "<<Tminmax[0]<<" "<<Tminmax[1]<<" "<<Iminmax[0]<<" "<<Iminmax[1]<<" "<<ICminmax[0]<<" "<<ICminmax[1]<<std::endl;
  std::cout<<" mean sigma "<<Imeansigma[0]<<"+-"<<Imeansigma[1]<<" "<<Energymeansigma[0]<<" "<<Energymeansigma[1]<<std::endl;
  
  std::cout<<" TestRange "<<RangeAll[1][1][0]<<" "<<RangeAll[1][3][0]<<std::endl;
  std::cout<<" type : "<<type<<" "<<ICtype<<std::endl;
  
  // TH3F* h_hist = new TH3F("h_hist","h_hist",NbinX,minHist,1.,NbinX,minHist,1.,NbinY,minHist,1.);
  // TH3F* h_histNameBeam = new TH3F("h_histNameBeam","h_histNameBeam",NbinX,minHist,1.,NbinX,minHist,1.,NbinY,minHist,1.);
  // TH3F* h_histCmTarget = new TH3F("h_histTickness","h_histThickness",NbinX,minHist,1.,NbinX,minHist,1.,NbinY,minHist,1.);
  // TH3F* h_histNameTarget = new TH3F("h_histNameTarget","h_histNameTarget",NbinX,minHist,1.,NbinX,minHist,1.,NbinY,minHist,1.);

  Int_t NbinX = Nbin;
  Int_t NbinY = Nbin;
  TH2F* h_hist = new TH2F("h_hist","h_hist",NbinX,minHist,1.,NbinY,minHist,1.);
  TH2F* h_histNameBeamPrimary = new TH2F("h_histNameBeamPrimary","h_histNameBeamPrimary",NbinX,minHist,1.,NbinY,minHist,1.);
  TH2F* h_histNameBeamSecondary = new TH2F("h_histNameBeamSecondary","h_histNameBeamSecondary",NbinX,minHist,1.,NbinY,minHist,1.);
  TH2F* h_histCrossSection = new TH2F("h_histCrosssection","h_histCrossSection",NbinX,minHist,1.,NbinY,minHist,1.);
  TH2F* h_histEnergy = new TH2F("h_histEnergy","h_histEnergy",NbinX,minHist,1.,NbinY,minHist,1.);
  TH2F* h_histCmTarget = new TH2F("h_histThickness","h_histThickness",NbinX,minHist,1.,NbinY,minHist,1.);
  TH2F* h_histNameTarget = new TH2F("h_histNameTarget","h_histNameTarget",NbinX,minHist,1.,NbinY,minHist,1.);

  TH2F* h_Internal = new TH2F("h_internal","h_internal",15,0,15,1000,-10,10);

  double maxMax = -99999.;
  std::set<dataMax,CompdataMax> maxMaxId;
  //std::vector<double> maxMaxIdValue(2);

  for(int i=1;i<=h_hist->GetNbinsX();++i)
    {
      for(int j=1;j<=h_hist->GetNbinsY();++j)
	//for(int k=1;k<=h_hist->GetNbinsZ();++k)
	{
	  double alpha = h_hist->GetXaxis()->GetBinCenter(i);
	  double beta  = h_hist->GetYaxis()->GetBinCenter(j);
	  //double gamma  = h_hist->GetZaxis()->GetBinCenter(k);

	  //double gamma = ComputeWeigth(alpha,beta,NtypeComp);

	  if(compute.ParAcceptable(alpha,beta,Ntype) == false)
	    continue;

	  double max_G = -1.e100;
	  int max_CmTarget = -1;
	  std::vector<int> max_nameIds(6,0);
	  double maxIntensity = 0.;
	  double maxContamination = 0.;
	  double maxEnergyM = 0.;
	  double maxEnergyS = 0.;
	  double maxSurvival = 0.;
	  double maxCrossSection = 0;

	  for(auto it_secondaryReaction = hyp_prod.cbegin(), it_secondaryReactionEnd = hyp_prod.cend(); it_secondaryReaction != it_secondaryReactionEnd; ++it_secondaryReaction)
	    {
	      
	      auto it_findStable = StableBeam.find(it_secondaryReaction->nameB);
	      if(it_findStable != StableBeam.end())
		continue;
	      auto it_findTargetS = SecondTargetDensity.find(it_secondaryReaction->nameT);
	      if(it_findTargetS == SecondTargetDensity.end())
		continue;
	      double SecondTargetD = it_findTargetS->second;

	      //std::cout<<"case:"<<it_secondaryReaction->nameB<<" "<<it_secondaryReaction->nameT<<" CX :["<<it_secondaryReaction->cross_section<<" ]"<<std::endl;
	      for(auto it_FBT = BeamTarget.cbegin(), it_FBT_end = BeamTarget.cend();it_FBT!= it_FBT_end;++it_FBT)
		{
		  auto it_findRejected = RejectedBeam.find(it_FBT->second.F);
		  if(it_findRejected != RejectedBeam.end())
		    continue;
		  //std::cout<<" -primary:"<<it_FBT->second.F<<" "<<it_FBT->first<<std::endl;
		  if(it_secondaryReaction->nameB == it_FBT->second.F )
		    {
		      std::vector<int> tempIds(it_FBT->second.BT_ids);
		    
		      for(unsigned int idCmTarget = 0; idCmTarget< it_FBT->second.TargetCm.size();++idCmTarget)
			{
			  //std::cout<<" idCm :"<<idCmTarget<<std::endl;
			  double tempT = it_FBT->second.TargetGramCm[idCmTarget];
			  double tempI = it_FBT->second.ProdFrag[idCmTarget];
			  if(tempI>Iminmax[1])
			    tempI=Iminmax[1];
			  double tempEpar[2] = {it_FBT->second.EnergyMean[idCmTarget],it_FBT->second.EnergySigma[idCmTarget]}; 

			  double tempICX = tempI*4.*SecondTargetD*0.2*it_secondaryReaction->cross_section*6.02214129e-07*StrangenessProdNorm(tempEpar[0]); // 6.02409638554217e-07 = 1e-30/1.66e-24
			  double tempC = it_FBT->second.ProdPara1[idCmTarget]+it_FBT->second.ProdPara2[idCmTarget]+it_FBT->second.ProdPara3[idCmTarget];
			  //double tempE = (2.e3-tempEpar[0])*(2.e3-tempEpar[0])/tempEpar[1]/tempEpar[1];
			  if(ICtype == "I-C")
			    tempC = tempI - tempC;
			  // if(it_FBT->second.BT_ids[0] == 14 && it_FBT->second.BT_ids[1] == 7 && it_FBT->second.BT_ids[2] == 9 && it_FBT->second.BT_ids[3] == 4 && it_FBT->second.BT_ids[4] == 12 && it_FBT->second.BT_ids[5] == 7 && TMath::Abs(it_FBT->second.TargetCm[idCmTarget]-5.)<1.e-3 && it_secondaryReaction->nameT=="C12")
			  //   {
			  //     std::cout<<"Test: "<<it_FBT->first<<std::endl;
			  //     it_FBT->second.Print(5.);
			  //     std::cout<<" Secondary: "<<it_secondaryReaction->nameB<<" "<<it_secondaryReaction->nameT<<" "<<tempICX<<std::endl;
			  //   }

			  //h_Internal->Fill(5.,tempICX/640.);

			  if(tempICX<0.5)
			    continue;

			  double Brho = compute.Brho(it_FBT->second.BT_ids[5],it_FBT->second.BT_ids[4],tempEpar[0]); 
			  h_Internal->Fill(14.,Brho/10.);
 
			  if(Brho>20.5)
			    continue;
			  //double tempE = 10*TMath::Gaus(2000,tempEpar[0],tempEpar[1]*5);
			  
			  //h_Internal->Fill(11.,(tempT-Tmeansigma[0])/(Tmeansigma[1]));
			  //h_Internal->Fill(12.,(tempI-Imeansigma[0])/(Imeansigma[1]));
			  //h_Internal->Fill(13.,(tempEpar[0]-Energymeansigma[0])/(Energymeansigma[1]));
			  
			  
			  std::vector<double> Input = {tempT,tempI,tempC,tempEpar[0]};
			  //std::vector<double> Input2 = {0,tempICX};
			  compute.Invariant(Input, RangeAll[type], type, h_Internal);
			  //Invariant(Input2, RangeAll[type], type);
			  
			  // h_Internal->Fill(6.,(tempT-Tminmax[0])/(Tminmax[1]-Tminmax[0]));
			  // h_Internal->Fill(7.,(tempI-Iminmax[0])/(Iminmax[1]-Iminmax[0]));
			  // if(ICtype == "I-C")
			  //   h_Internal->Fill(8.,(tempC-ICminma[0])/(ICminmax[1]-ICminmax[0]));
			  // else
			  //   h_Internal->Fill(8.,(tempC-Cminmax[0])/(Cminmax[1]-Cminmax[0]));
			
			
			  // if(type==0)
			  //   {
			  //     tempT -= Tminmax[0];
			  //     tempT /= Tminmax[1]-Tminmax[0];
			    
			  //     tempI -= Iminmax[0];
			  //     tempI /= Iminmax[1]-Iminmax[0];
			    
			  //     if(ICtype == "I-C")
			  //       {
			  // 	tempC -= ICminmax[0];
			  // 	tempC /= ICminmax[1]-ICminmax[0];
			  //       }
			  //     else
			  //       {
			  //       tempC -= Cminmax[0];
			  //       tempC /= Cminmax[1]-Cminmax[0];
			  //       }
			  //   }
			  // else
			  //   {
			  //     tempT -= Tmeansigma[0];
			  //     tempT /= Tmeansigma[1];
			    
			  //     tempI -= Imeansigma[0];
			  //     tempI /= Imeansigma[1];
			    
			  //     if(ICtype=="I-C")
			  //       {
			  // 	tempC -= ICmeansigma[0];
			  // 	tempC /= ICmeansigma[1];
			  //       }
			  //     else
			  //       {
			  // 	tempC -= Cmeansigma[0];
			  // 	tempC /= Cmeansigma[1];
			  //       }
			  //   }
			  // double norm = TMath::Abs(alpha + beta - gamma);
			  // double tempG = alpha*tempI + beta*tempT - gamma*tempC ;
			  // tempG /= norm;
			
			  //double tempG = alpha*tempI - beta*tempT + (1-alpha+beta)*tempC ;
			  std::vector<double> inPar(2);
			  inPar[0] = alpha;
			  inPar[1] = beta;

			  std::vector<double> outPar(3);
			  int res = compute.setPermutation(Mtype,inPar,outPar,NtypeComp);
			  if(res==-1)
			    return -1;
			  if(res==-2)
			    continue;
			  //double tempG =  - 1./0.5808*1./0.54*1/1.10*outPar[0]*Input[0] + 1./320*1./0.32*outPar[1]*tempICX + outPar[2]*Input[3] ;//+ outPar[2]*Input[2] ;
			  double tempG =  - 1./0.5808*1./0.54*1/1.10*1./0.92*outPar[0]*Input[0] + (1./320*1./0.32*1./0.14*1./1.02*1./0.98)*outPar[1]*tempICX + 1./0.98*outPar[2]*((Input[3]-0.46)/(1.-0.46)) + 1./1.02*0.5*Input[1];//+ outPar[2]*Input[2]
			  //double tempG =  - 1./0.5808*outPar[0]*Input[0] + 1./0.56448*outPar[1]*Input[1] + outPar[2]*Input[3] + 1./5.4*Input2[1];
		  
			  // h_Internal->Fill(0.,tempT);
			  // h_Internal->Fill(1.,tempI);
			  // h_Internal->Fill(2.,tempC);
			  // h_Internal->Fill(3.,tempEpar[0]);
			  // h_Internal->Fill(4.,tempG);
			  h_Internal->Fill(0.,1./0.5808*1./0.54*1./1.10*1./0.92*Input[0]);
			  h_Internal->Fill(1.,1./1.02*Input[1]);
			  h_Internal->Fill(2.,1./320.*1./0.32*1./0.14*1./1.02*1./0.98*tempICX);
			  h_Internal->Fill(3.,1./0.98*(Input[3]-0.46)/(1.-0.46));

			  if(tempG > max_G)
			    {
			      max_G = tempG;
			      max_CmTarget = it_FBT->second.TargetCm[idCmTarget];//idCmTarget;
			      max_nameIds = tempIds;
			      maxIntensity = it_FBT->second.ProdFrag[idCmTarget];
			      maxContamination = it_FBT->second.ProdPara1[idCmTarget]+it_FBT->second.ProdPara2[idCmTarget]+it_FBT->second.ProdPara3[idCmTarget];
			      maxEnergyM = it_FBT->second.EnergyMean[idCmTarget];
			      maxEnergyS = it_FBT->second.EnergySigma[idCmTarget];
			      maxSurvival = it_FBT->second.Survival[idCmTarget];
			      maxCrossSection = tempICX;//it_FBT->second.ProdFrag[idCmTarget]*it_secondaryReaction->cross_section;
			    }
			}
		    }
		}
	    }

	  double nameTemplateB = max_nameIds[0]*100 + max_nameIds[1];
	  double nameTemplateT = max_nameIds[2]*100 + max_nameIds[3];
	  double nameTemplateF = max_nameIds[4]*100 + max_nameIds[5];
	  
	  // h_hist->Fill(alpha,beta,gamma,max_G);
	  // h_histNameBeam->Fill(alpha,beta,gamma,nameTemplateT);
	  // h_histCmTarget->Fill(alpha,beta,gamma,max_CmTarget);
	  // h_histNameTarget->Fill(alpha,beta,gamma,nameTemplateB);
	  if(TMath::Abs(max_G-maxMax)<=1e-4)
	    {
	      
	      dataMax tempMax = {i,j,alpha,beta, nameTemplateB, nameTemplateT, nameTemplateF, max_CmTarget, maxIntensity, maxContamination, maxEnergyM, maxEnergyS,maxSurvival, maxCrossSection};
	      //maxMaxId.push_back(tempMax);
	      maxMaxId.insert(tempMax);
	      
	      //std::cout<<" max~ #"<<i<<" "<<j<<" "<<max_G<<" | "<<maxMax <<" size:"<<maxMaxId.size()<<std::endl;
	    }
	  else if(max_G>maxMax)
	    {
	      dataMax tempMax = {i,j, alpha, beta, nameTemplateB, nameTemplateT, nameTemplateF, max_CmTarget, maxIntensity, maxContamination, maxEnergyM, maxEnergyS,maxSurvival, maxCrossSection};
	      //maxMaxId[0] = i;
	      //maxMaxId[1] = j;

	      //maxMaxIdValue[0] = alpha;
	      //maxMaxIdValue[1] = beta;
	      if(maxMaxId.size()==0)
		//maxMaxId.push_back(tempMax);
		maxMaxId.insert(tempMax);
	      else
		{
		  maxMaxId.clear();
		  maxMaxId.insert(tempMax);
		  //maxMaxId.push_back(tempMax);
		}
	      //std::cout<<" Max #"<<i<<" "<<j<<" "<<max_G<<" | "<<maxMax<<" size:"<<maxMaxId.size()<<std::endl;
	      
	      maxMax = max_G;
	    }

	  
	  h_hist->Fill(alpha,beta,max_G);
	  h_histNameBeamPrimary->Fill(alpha,beta,nameTemplateB);
	  h_histNameBeamSecondary->Fill(alpha,beta,nameTemplateF);
	  h_histCmTarget->Fill(alpha,beta,max_CmTarget);
	  h_histNameTarget->Fill(alpha,beta,nameTemplateT);
	  h_histEnergy->Fill(alpha,beta,TMath::FloorNint(maxEnergyM*10)*1e-4);
	  h_histCrossSection->Fill(alpha,beta,TMath::FloorNint(maxCrossSection*10)*0.1);


	}
    }
  //std::cout<<" Max of the maxG : ["<<maxMaxId[0]<<", "<<maxMaxId[1]<<"] | alpha="<<maxMaxIdValue[0]<<" beta="<<maxMaxIdValue[1]<<std::endl;
  std::cout<<" #max :"<<maxMaxId.size()<<std::endl;
  std::vector<TMarker*> maxMarkers (maxMaxId.size(),0);

  //for(unsigned int maxIndex = 0; maxIndex < maxMaxId.size(); ++maxIndex)
  int maxIndex = 0;
  for(auto dataId : maxMaxId)
    {
      std::cout<<" Max ["<<maxIndex<<"] -> ("<<dataId.idI<<", "<<dataId.idJ<<") | alpha="<<dataId.alpha<<" beta="<<dataId.beta<<" "<<std::endl;
      std::cout<<" +-> beamP:"<<dataId.beam<<" target:"<<dataId.target<<" IDcm:"<<dataId.cmtarget<<" Int:"<<dataId.intensity<<" | cont:"<<dataId.contamination<<std::endl;
      std::cout<<" +-> beamS:"<<dataId.fragment<<" Energy:"<<dataId.energyMean<<" +- "<<dataId.energySigma<<" Survival:"<<dataId.survivalRate<<std::endl;
      std::cout<<" +-> CX:"<<dataId.crossSectionHyp<<std::endl;
      maxMarkers[maxIndex] = new TMarker(dataId.alpha,dataId.beta,30);
      ++maxIndex;
    }

  std::string OptionDraw = NbinX > 30 || NbinY > 30 ? "colz" : "colz text";
  
  TCanvas* c1 = new TCanvas("Internal","Internal",500,500);
  c1->cd();
  h_Internal->Draw("colz");

  TCanvas* c2 = new TCanvas("Cost","Cost",500,500);
  //c2->Divide(2,2);
  // c2->cd(1);
  // h_hist->Project3D("xy")->Draw("colz");
  // c2->cd(2);
  // h_hist->Project3D("xz")->Draw("colz");
  // c2->cd(3);
  // h_hist->Project3D("yz")->Draw("colz");
  //c2->cd(4);
  h_hist->Draw("colz");
  for(unsigned int maxIndex = 0; maxIndex < maxMarkers.size(); ++maxIndex)
    maxMarkers[maxIndex]->Draw("same");
  

  TCanvas* c3 = new TCanvas("NameBeam","NameBeam",500,500);
  //c3->Divide(2,2);
  // c3->cd(1);
  // h_histNameBeam->Project3D("xy")->Draw("colz");  
  // c3->cd(2);
  // h_histNameBeam->Project3D("xz")->Draw("colz");
  // c3->cd(3);
  // h_histNameBeam->Project3D("yz")->Draw("colz");
  //c3->cd(4);
  h_histNameBeamPrimary->Draw(OptionDraw.c_str());
  for(unsigned int maxIndex = 0; maxIndex < maxMarkers.size(); ++maxIndex)
    maxMarkers[maxIndex]->Draw("same");

  TCanvas* c3_1 = new TCanvas("NameTarget","NameTarget",500,500);
  //c3_1->Divide(2,2);
  // c3_1->cd(1);
  // h_histNameTarget->Project3D("xy")->Draw("colz");
  // c3_1->cd(2);
  // h_histNameTarget->Project3D("xz")->Draw("colz");
  // c3_1->cd(3);
  // h_histNameTarget->Project3D("yz")->Draw("colz");
  //c3_1->cd(4);
  h_histNameTarget->Draw(OptionDraw.c_str());
  for(unsigned int maxIndex = 0; maxIndex < maxMarkers.size(); ++maxIndex)
    maxMarkers[maxIndex]->Draw("same");

  TCanvas* c3_2 = new TCanvas("NameBeam2","NameBeam2",500,500);
  //c3->Divide(2,2);
  // c3->cd(1);
  // h_histNameBeam->Project3D("xy")->Draw("colz");  
  // c3->cd(2);
  // h_histNameBeam->Project3D("xz")->Draw("colz");
  // c3->cd(3);
  // h_histNameBeam->Project3D("yz")->Draw("colz");
  //c3->cd(4);
  h_histNameBeamSecondary->Draw(OptionDraw.c_str());
  for(unsigned int maxIndex = 0; maxIndex < maxMarkers.size(); ++maxIndex)
    maxMarkers[maxIndex]->Draw("same");

  TCanvas* c4 = new TCanvas("Thickness","Thickness",500,500);
  //c4->Divide(2,2);
  //c4->cd(1);
  // h_histCmTarget->Project3D("xy")->Draw("colz");
  // c4->cd(2);
  // h_histCmTarget->Project3D("xz")->Draw("colz");
  // c4->cd(3);
  // h_histCmTarget->Project3D("yz")->Draw("colz");
  //c4->cd(4);
  h_histCmTarget->Draw(OptionDraw.c_str());
  for(unsigned int maxIndex = 0; maxIndex < maxMarkers.size(); ++maxIndex)
    maxMarkers[maxIndex]->Draw("same");

  TCanvas* c5 = new TCanvas("CX","CX",500,500);
  //c3->Divide(2,2);
  // c3->cd(1);
  // h_histNameBeam->Project3D("xy")->Draw("colz");  
  // c3->cd(2);
  // h_histNameBeam->Project3D("xz")->Draw("colz");
  // c3->cd(3);
  // h_histNameBeam->Project3D("yz")->Draw("colz");
  //c3->cd(4);
  h_histCrossSection->Draw(OptionDraw.c_str());
  for(unsigned int maxIndex = 0; maxIndex < maxMarkers.size(); ++maxIndex)
    maxMarkers[maxIndex]->Draw("same");

  TCanvas* c6 = new TCanvas("Energy","Energy",500,500);
  //c3->Divide(2,2);
  // c3->cd(1);
  // h_histNameBeam->Project3D("xy")->Draw("colz");  
  // c3->cd(2);
  // h_histNameBeam->Project3D("xz")->Draw("colz");
  // c3->cd(3);
  // h_histNameBeam->Project3D("yz")->Draw("colz");
  //c3->cd(4);
  h_histEnergy->Draw(OptionDraw.c_str());
  for(unsigned int maxIndex = 0; maxIndex < maxMarkers.size(); ++maxIndex)
    maxMarkers[maxIndex]->Draw("same");



  return 0;
}


int AnalysisHypStable(std::string namefileEpax, std::string namefileHyp, std::string Hyp)
{
  
  //"n","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr",
  // 0 , 1 , 2  , 3  , 4  , 5 , 6 , 7 , 8 , 9 , 10 , 11 , 12 , 13 , 14 , 15, 16, 17 , 18 , 19, 20 , 21 , 22 , 23, 24 ,

  const std::map<int,int> StableBeam = {{1,1},{1,2},{2,3},{2,4},{3,6},{3,7},{4,9},{5,10},{5,11},{6,12},{6,13},{7,14},{7,15},{8,16},{8,17},{8,18},{9,19},{10,20},{10,21},{10,22},{11,23},{12,24},{12,25},{12,26},{13,27},{14,28},{14,29},{14,30},{15,31},{16,32},{16,33},{16,34},{16,36},{17,35},{17,37},{18,36},{18,38},{18,40}};


  
  std::map<std::string,double> SecondTargetDensity;
  SecondTargetDensity["C12"] = 2.21/12;
  //SecondTargetDensity["Be9"] = 1.85/9;
    
  double PrimaryBeamInt = 1e10;

  std::ifstream ifs ( namefileEpax.c_str() );

  //std::string temp_line;
  int Af,Zf,Ap,Zp,At,Zt;
  double TargetDensity,xs,ProdRate;
  int Sister1,Sister2,Sister3,ZZ;
  double CmTarget, gCmTarget;
  double SurvivalRate, MeanEnergy, SigmaEnergy, Transmission;
  double prod1;
  double prod2;
  double prod3;
  double prod4;
  
  std::map<std::string,tuple_graph > BeamTarget;

  std::set<std::string> FragAvailable;
  
  std::string temp_line;
  std::getline(ifs,temp_line);
  
  while(std::getline(ifs,temp_line))
    {
      std::stringstream stream(temp_line);
      stream >> Af >> Zf >> Ap >> Zp >> At >> Zt >> TargetDensity >> xs >> ProdRate >> Sister1 >> Sister2 >> Sister3 >> ZZ
	     >> CmTarget >> SurvivalRate >> MeanEnergy >> SigmaEnergy >> Transmission >> prod1 >> prod2 >> prod3 >> prod4 ; 

      if(MeanEnergy<700.)
	continue;
      if(At == 2 && Zt == 1)
	continue;
      if(MeanEnergy-3*SigmaEnergy < 1600.)
	continue;

      // auto it_Z = StableBeam.find(Zf);
      // if(it_Z == StableBeam.end())
      // 	continue;
      // if(Af != it_Z->second+1 && Af != it_Z->second-1)
      // 	continue;
      
      std::vector<int> idsBT ;
      idsBT.push_back(Ap);
      idsBT.push_back(Zp);
      idsBT.push_back(At);
      idsBT.push_back(Zt);
      idsBT.push_back(Af);
      idsBT.push_back(Zf);
      TString nameF("");
      nameF+= ElName2[Zf];
      nameF+=Af;
      TString nameP("");
      nameP+= ElName2[Zp];
      nameP+=Ap;
      TString nameT(""); 
      nameT+= ElName2[Zt];
      nameT+=At;

      std::string SnameF(nameF.Data());
      std::string SnameP(nameP.Data());
      std::string SnameT(nameT.Data());
      std::string SnameAll(SnameF);
      SnameAll+=SnameP;
      SnameAll+=SnameT;
      std::string SnameBT(SnameP);
      SnameBT+="+";
      SnameBT+=SnameT;

      FragAvailable.emplace(SnameF);
      
      std::map<std::string,tuple_graph>::iterator it_FBT = BeamTarget.find(SnameAll);
      if(it_FBT==BeamTarget.end())
	{
	  tuple_graph temp(SnameBT,SnameF,idsBT,CmTarget,ProdRate*CmTarget,PrimaryBeamInt*prod1*(1.-SurvivalRate),PrimaryBeamInt*prod2*(1.-SurvivalRate),PrimaryBeamInt*prod3*(1.-SurvivalRate),PrimaryBeamInt*prod4*(1.-SurvivalRate),MeanEnergy,SigmaEnergy,SurvivalRate);
	  BeamTarget.insert(std::pair<std::string,tuple_graph>(SnameAll,temp));
	}
      else
	{
	  it_FBT->second.AddValue(CmTarget,ProdRate*CmTarget,PrimaryBeamInt*prod1*(1.-SurvivalRate),PrimaryBeamInt*prod2*(1.-SurvivalRate),PrimaryBeamInt*prod3*(1.-SurvivalRate),PrimaryBeamInt*prod4*(1.-SurvivalRate),MeanEnergy,SigmaEnergy,SurvivalRate);
	}      
    }

  for(auto& it_stable : StableBeam)
    {
      
      std::vector<int> idsBT ;
      idsBT.push_back(0);
      idsBT.push_back(0);
      idsBT.push_back(0);
      idsBT.push_back(0);
      idsBT.push_back(it_stable.second);
      idsBT.push_back(it_stable.first);
      TString nameF("");
      nameF+= ElName2[it_stable.first];
      nameF+=it_stable.second;
      TString nameP("Primary");
      TString nameT("Primary"); 

      std::string SnameF(nameF.Data());
      std::string SnameP(nameP.Data());
      std::string SnameT(nameT.Data());
      std::string SnameAll(SnameF);
      SnameAll+=SnameP;
      SnameAll+=SnameT;
      std::string SnameBT(SnameP);
      SnameBT+="+";
      SnameBT+=SnameT;

      FragAvailable.emplace(SnameF);

      std::map<std::string,tuple_graph>::iterator it_FBT = BeamTarget.find(SnameAll);
      if(it_FBT==BeamTarget.end())
	{
	  tuple_graph temp(SnameBT,SnameF,idsBT,0,0,PrimaryBeamInt,0,0,0,2000.,0.1,1);
	  BeamTarget.insert(std::pair<std::string,tuple_graph>(SnameAll,temp));
	}
      else
	{
	  std::cout<<"!> stable beam should be unique ! "<<SnameAll<<std::endl;
	}      
    }
  
  
  std::ifstream ifsHyp (namefileHyp.c_str());
  std::set<nameBTtuple,CompBTTuple> hyp_prodAll;

  std::getline(ifsHyp,temp_line);
  double max_CX = -1 ;
  while(std::getline(ifsHyp,temp_line))
    {
      std::stringstream stream(temp_line);
      std::string hyp;
      int nb;
      stream >> hyp >> nb;
      if(hyp != Hyp)
	continue;
      std::cout<<temp_line<<std::endl;
      nb /= 2;
      for(int id = 0; id<nb ; ++id)
	{
	  std::string BT;
	  double CX1;

	  stream >> BT >> CX1;

	  std::cout<<" ---> id#"<<id<<" "<<BT<<" "<<CX1;
	  std::size_t foundPlus = BT.find("+");
	  if(foundPlus == std::string::npos )
	    {
	      std::cout<<"E> no + ! "<<id<<" "<<BT<<" "<<hyp<<std::endl;
	      return -1;
	    }
	  nameBTtuple nameBT;
	  nameBT.nameB=BT.substr(0,foundPlus);
	  nameBT.nameT=BT.substr(foundPlus+1);
	  nameBT.cross_section = CX1;

	  if(nameBT.nameT!="C12")
	    continue;
	  // auto isSecondaryBeamAvailable = FragAvailable.find(nameBT.nameB);
	  // if(isSecondaryBeamAvailable != FragAvailable.end())
	  //   {
	  std::cout<<" beam:"<<nameBT.nameB<<" target:"<<nameBT.nameT<<" "<<nameBT.cross_section<<std::endl; 
	  auto ret = hyp_prodAll.emplace(nameBT);
	  if (!ret.second)
	    std::cout << "already exists in hyp_prod\n";
	  else
	    {
	      if(CX1 > max_CX )
		max_CX = CX1;
	    }
	  //   }
	  // else
	  //   std::cout<<" 2nd beam not available !"<<" beam:"<<nameBT.nameB<<" target:"<<nameBT.nameT<<" "<<nameBT.cross_section<<std::endl;
	}
    }

  auto Pred = [&] (const nameBTtuple& nameBT) -> bool
	       {
		 auto isFound = FragAvailable.find(nameBT.nameB);
		 if(isFound != FragAvailable.end())
		   return true;
		 else
		   {
		     std::cout<<" 2nd beam not available !"<<" beam:"<<nameBT.nameB<<" target:"<<nameBT.nameT<<" "<<nameBT.cross_section<<" "<<nameBT.cross_section/max_CX<<std::endl;
		     return false;
		   }
	       };
  
  std::set<nameBTtuple,CompBTTuple> hyp_prod;
  for(auto item : hyp_prodAll)
    if(Pred(item))
      {
	hyp_prod.emplace(item);
      }
  
  
  if(hyp_prod.size()==0)
    {
      std::cout<<"!> no hypernuclei selected ! "<<Hyp<<std::endl;
      return -2;
    }

  std::cout<<"Hyp required: "<<Hyp<<std::endl<<" 2nd Reaction:";

  // for(std::map<std::string, tuple_graph>::const_iterator it_FBT = BeamTarget.begin(), it_FBT_end = BeamTarget.end();it_FBT!= it_FBT_end;++it_FBT)
  //   {
  //      std::cout<<"*Entry :"<<it_FBT->first<<std::endl;
  //      it_FBT->second.Print();
  //   }

  for(auto& FragA : FragAvailable)
    {
      std::cout<<FragA<<std::endl;
    }
  
  for(auto& hypBT : hyp_prod)
    {
      //if(hypBT.nameT=="C12" && hypBT.nameB=="N12")
      std::cout<<" <"<<hypBT.nameB<<"> + <"<<hypBT.nameT<<"> ("<<hypBT.cross_section<<")"<<std::endl;
    }
  std::cout<<std::endl;


  for(std::map<std::string, tuple_graph>::const_iterator it_FBT = BeamTarget.begin(), it_FBT_end = BeamTarget.end();it_FBT!= it_FBT_end;++it_FBT)
    {
      if(it_FBT->second.F == "C11") 
	{
	  std::cout<<"*Entry :"<<it_FBT->first<<std::endl;
	  if(/*it_FBT->second.BT_ids[0] == 14 && it_FBT->second.BT_ids[1] == 7 &&*/ it_FBT->second.BT_ids[2] == 9 && it_FBT->second.BT_ids[3] == 4)
	    it_FBT->second.Print();
	}
    }
  
  

}
