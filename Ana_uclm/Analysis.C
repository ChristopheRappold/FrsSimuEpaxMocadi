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


auto AnalysisHyp(std::string namefileEpax, std::string namefileHyp, std::string Hyp, int type = 0, std::string ICtype = "I-C", std::string Mtype = "A B AB", int Ntype = -1, int NtypeComp = 1, double minHist = 0, int Nbin = 20, bool interp = false, bool withStableBeam = false, bool noHisto = false)
{
  std::cout<<" Processing :"<<Hyp<<" "<<noHisto<<std::endl;

  const std::set<std::string> StableBeam = {"H1","H2","He3","He4","Li6","Li7","Be9","B10","B11","C12","C13","N14","N15","O16","O17","O18","F19","Ne20","Ne21","Ne22","Na23","Mg24","Mg25","Mg26","Al27","Si28","Si29","Si30","P31","S32","S33","S34","S36","Cl35","Cl37","Ar36","Ar38","Ar40"};
  
  std::map<std::string,double> SecondTargetDensity;
  SecondTargetDensity["C12"] = 2.21/12;
  SecondTargetDensity["Be9"] = 1.85/9;
  
  const double PrimaryBeamInt = 5.e9;
  const double IntensityMax = 5.e6;
  
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
  
  std::map<std::string,tuple_graph > BeamTargetOriginal;

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
      
      std::map<std::string,tuple_graph>::iterator it_FBT = BeamTargetOriginal.find(SnameAll);
      if(it_FBT==BeamTargetOriginal.end())
	{
	  tuple_graph temp(SnameBT,SnameF,idsBT,CmTarget,ProdRate*CmTarget,PrimaryBeamInt*prod1*(1.-SurvivalRate),PrimaryBeamInt*prod2*(1.-SurvivalRate),PrimaryBeamInt*prod3*(1.-SurvivalRate),PrimaryBeamInt*prod4*(1.-SurvivalRate),MeanEnergy,SigmaEnergy,SurvivalRate,Transmission);
	  BeamTargetOriginal.insert(std::pair<std::string,tuple_graph>(SnameAll,temp));
	}
      else
	{
	  it_FBT->second.AddValue(CmTarget,ProdRate*CmTarget,PrimaryBeamInt*prod1*(1.-SurvivalRate),PrimaryBeamInt*prod2*(1.-SurvivalRate),PrimaryBeamInt*prod3*(1.-SurvivalRate),PrimaryBeamInt*prod4*(1.-SurvivalRate),MeanEnergy,SigmaEnergy,SurvivalRate,Transmission);
	}

      
      // if(IntCon.back() < MinMaxCheck[0])
      // 	{
      // 	  MinMaxCheck[0] = IntCon.back();
      // 	  NameCheck[0] = SnameAll;
      // 	}
      // if(IntCon.back() > MinMaxCheck[1])
      // 	{
      // 	  MinMaxCheck[1] = IntCon.back();
      // 	  NameCheck[1] = SnameAll;
      // 	}
      
    }



  
  std::ifstream ifsHyp (namefileHyp.c_str());
  std::set<HypDataProd,CompHypDataProd> hyp_prodAll;

  //std::getline(ifsHyp,temp_line);
  double max_CX = -1 ;
  while(std::getline(ifsHyp,temp_line))
    {
      std::stringstream stream(temp_line);
      std::string hyp;
      int nb;
      stream >> hyp >> nb;
      if(hyp != Hyp)
	continue;
      //std::cout<<temp_line<<std::endl;
      nb /= 2;
      for(int id = 0; id<nb ; ++id)
	{
	  std::string BT;
	  double CX1;

	  stream >> BT >> CX1;

	  //std::cout<<" ---> id#"<<id<<" "<<BT<<" "<<CX1;
	  std::size_t foundPlus = BT.find("+");
	  if(foundPlus == std::string::npos )
	    {
	      std::cout<<"E> no + ! "<<id<<" "<<BT<<" "<<hyp<<std::endl;
	      return std::make_tuple(static_cast<size_t>(0),0.,0.,-1.,-1.,0.,-1.,-1.);
	    }
	  HypDataProd nameBT;
	  nameBT.nameB=BT.substr(0,foundPlus);
	  nameBT.nameT=BT.substr(foundPlus+1);
	  nameBT.cross_section = CX1;

	  auto it_findTargetS = SecondTargetDensity.find(nameBT.nameT);
	  if(it_findTargetS == SecondTargetDensity.end())
	    continue;
	  nameBT.SecondTargetD = it_findTargetS->second;

	  // auto isSecondaryBeamAvailable = FragAvailable.find(nameBT.nameB);
	  // if(isSecondaryBeamAvailable != FragAvailable.end())
	  //   {
	  //std::cout<<" beam:"<<nameBT.nameB<<" target:"<<nameBT.nameT<<" "<<nameBT.cross_section<<std::endl; 
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

    
  std::map<std::string,tuple_graph > BeamTarget;
  if(interp == true)
    {
      Interpole doInterpoling;
      for(auto it_FBT = BeamTargetOriginal.cbegin(), it_FBT_end = BeamTargetOriginal.cend(); it_FBT != it_FBT_end; ++it_FBT)
	{
	  auto tempTuple = doInterpoling.Inter(it_FBT->second);
	  BeamTarget.emplace(it_FBT->first,tempTuple);
	}
    }
  else
    BeamTarget = BeamTargetOriginal;

  //const std::map<int,int> StableBeam = {{1,1},{1,2},{2,3},{2,4},{3,6},{3,7},{4,9},{5,10},{5,11},{6,12},{6,13},{7,14},{7,15},{8,16},{8,17},{8,18},{9,19},{10,20},{10,21},{10,22},{11,23},{12,24},{12,25},{12,26},{13,27},{14,28},{14,29},{14,30},{15,31},{16,32},{16,33},{16,34},{16,36},{17,35},{17,37},{18,36},{18,38},{18,40}};
  if(withStableBeam==true)
    {
      const double BeamIntensityPrimary = 1.e7;
      for(auto& it_stable : std::map<int,int>({{1,1},{1,2},{2,3},{2,4},{3,6},{3,7},{4,9},{5,10},{5,11},{6,12},{6,13},{7,14},{7,15},{8,16},{8,17},{8,18},{9,19},{10,20},{10,21},{10,22},{11,23},{12,24},{12,25},{12,26},{13,27},{14,28},{14,29},{14,30},{15,31},{16,32},{16,33},{16,34},{16,36},{17,35},{17,37},{18,36},{18,38},{18,40}}) )
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
	      
	      tuple_graph temp(SnameBT,SnameF,idsBT,0,0,BeamIntensityPrimary,0,0,0,2000.,0.1,1.,1.);
	      BeamTarget.insert(std::pair<std::string,tuple_graph>(SnameAll,temp));
	    }
	  else
	    {
	      std::cout<<"!> stable beam should be unique ! "<<SnameAll<<std::endl;
	    }      
	}
    }
  // std::cout<<"FragAvailable :";
  // for(auto& item : FragAvailable)
  //   {
  //     std::cout<<" "<<item;
  //   }
  // std::cout<<std::endl;
  
  auto Pred = [&FragAvailable,noHisto,max_CX] (const HypDataProd& nameBT) -> bool
	       {
		 auto isFound = FragAvailable.find(nameBT.nameB);
		 if(isFound != FragAvailable.end())
		   return true;
		 else
		   {
		     if(noHisto==false)
		       std::cout<<" 2nd beam not available !"<<" beam:"<<nameBT.nameB<<" target:"<<nameBT.nameT<<" "<<nameBT.cross_section<<" "<<nameBT.cross_section/max_CX<<std::endl;
		     return false;
		   }
	       };
  
  std::set<HypDataProd,CompHypDataProd> hyp_prod;
  for(auto item : hyp_prodAll)
    if(Pred(item))
      {
	hyp_prod.emplace(item);
      }
  
  
  if(hyp_prod.size()==0)
    {
      std::cout<<"!> no hypernuclei selected ! "<<Hyp<<std::endl;
      return std::make_tuple(static_cast<size_t>(0),0.,0.,-1.,-1.,0.,-1.,-1.);
    }

  // std::cout<<"Hyp required: "<<Hyp<<std::endl<<" 2nd Reaction:";
  // for(auto hypBT : hyp_prod)
  //   {
  //     //if(hypBT.nameT=="C12" && hypBT.nameB=="N12")
  //     std::cout<<" <"<<hypBT.nameB<<"> + <"<<hypBT.nameT<<"> ("<<hypBT.cross_section<<")"<<std::endl;
  //   }
  // std::cout<<std::endl;

  

  // for(std::map<std::string, tuple_graph>::const_iterator it_FBT = BeamTarget.begin(), it_FBT_end = BeamTarget.end();it_FBT!= it_FBT_end;++it_FBT)
  //    {
  //      std::cout<<"*Entry :"<<it_FBT->first<<std::endl;
  //      it_FBT->second.Print();
  //    }

  // auto it_FBTmin = BeamTarget.find(NameCheck[0]);
  // std::cout<<"* Entry min :"<<it_FBTmin->first<<std::endl;
  // it_FBTmin->second.Print();

  // auto it_FBTmax = BeamTarget.find(NameCheck[1]);
  // std::cout<<"* Entry max :"<<it_FBTmax->first<<std::endl;
  // it_FBTmax->second.Print();

  if(noHisto==false)
    {
      auto it_FBTtest = BeamTarget.find("N12N14Be9");
      std::cout<<"* Entry test :"<<it_FBTtest->first<<std::endl;
      it_FBTtest->second.Print();
    }
  Compute compute;
  RejectedBeamStruct Reject;
  std::vector<double> intensity;
  std::vector<double> contamination;
  std::vector<double> targetLength;
  std::vector<double> EnergyM;
  std::vector<double> Survival;
  std::vector<double> Trans;
  std::vector<double> IntCon;
  std::vector<double> crossX;

  Minimizer minimizer(compute,ICtype, type, Mtype, NtypeComp, 2000.,0.5);
  
  for(auto it_secondaryReaction = hyp_prod.cbegin(), it_secondaryReactionEnd = hyp_prod.cend(); it_secondaryReaction != it_secondaryReactionEnd; ++it_secondaryReaction)
    {
      if(withStableBeam==false)
	{
	  auto it_findStable = StableBeam.find(it_secondaryReaction->nameB);
	  if(it_findStable != StableBeam.end())
	    continue;
	}
      auto it_findTargetS = SecondTargetDensity.find(it_secondaryReaction->nameT);
      if(it_findTargetS == SecondTargetDensity.end())
	continue;
      auto HypData = *it_secondaryReaction;
      //std::cout<<"case:"<<it_secondaryReaction->nameB<<" "<<it_secondaryReaction->nameT<<" CX :["<<it_secondaryReaction->cross_section<<" ]"<<std::endl;
      for(auto it_FBT = BeamTarget.cbegin(), it_FBT_end = BeamTarget.cend();it_FBT!= it_FBT_end;++it_FBT)
	{
	  auto it_findRejected = Reject.RejectedBeam.find(it_FBT->second.F);
	  if(it_findRejected != Reject.RejectedBeam.end() && it_FBT->second.BT !="Primary+Primary")
	    continue;

	  if(HypData.nameB == it_FBT->second.F )
	    {	
	      for(unsigned int idCmTarget = 0; idCmTarget< it_FBT->second.TargetCm.size();++idCmTarget)
		{

		  auto res_Min = minimizer.CheckConditions(0.4,0.6,idCmTarget,it_FBT->second,HypData,IntensityMax);
		  int error_res = std::get<0>(res_Min);
		  if(error_res!=0)
		    continue;
		  
		  intensity.push_back(it_FBT->second.ProdFrag[idCmTarget]);
		  contamination.push_back(it_FBT->second.ProdPara1[idCmTarget]+it_FBT->second.ProdPara2[idCmTarget]+it_FBT->second.ProdPara3[idCmTarget]);
		  targetLength.push_back(it_FBT->second.TargetGramCm[idCmTarget]);      
		  IntCon.push_back(it_FBT->second.ProdFrag[idCmTarget]-(it_FBT->second.ProdPara1[idCmTarget]+it_FBT->second.ProdPara2[idCmTarget]+it_FBT->second.ProdPara3[idCmTarget]));
		  EnergyM.push_back(it_FBT->second.EnergyMean[idCmTarget]);
		  Survival.push_back(1.-it_FBT->second.Survival[idCmTarget]);
		  Trans.push_back(it_FBT->second.Trans[idCmTarget]);
		  crossX.push_back(std::get<2>(res_Min));
		}
	    }
	}
    }
  if(crossX.size()==0)
    {
      std::cout<<"--- !> Nothing to search "<<std::endl;
      return std::make_tuple(static_cast<size_t>(0),0.,0.,-1.,-1.,0.,-1.,-1.);
    }
  Regularize Reg(targetLength, intensity, contamination, IntCon, EnergyM, Survival, Trans, crossX, ICtype, IntensityMax);

  if(noHisto == false)
    {
      std::cout<<" mix max "<<Reg.Tminmax[0]<<" "<<Reg.Tminmax[1]<<" "<<Reg.Iminmax[0]<<" "<<Reg.Iminmax[1]<<" "<<Reg.ICminmax[0]<<" "<<Reg.ICminmax[1]<<std::endl;
      std::cout<<" mean sigma "<<Reg.Imeansigma[0]<<"+-"<<Reg.Imeansigma[1]<<" "<<Reg.Energymeansigma[0]<<" "<<Reg.Energymeansigma[1]<<std::endl;
  
      std::cout<<" TestRange "<<Reg.RangeAll[1][1][0]<<" "<<Reg.RangeAll[1][3][0]<<std::endl;
      std::cout<<" type : "<<type<<" "<<ICtype<<std::endl;
    }

  // TH3F* h_hist = new TH3F("h_hist","h_hist",NbinX,minHist,1.,NbinX,minHist,1.,NbinY,minHist,1.);
  // TH3F* h_histNameBeam = new TH3F("h_histNameBeam","h_histNameBeam",NbinX,minHist,1.,NbinX,minHist,1.,NbinY,minHist,1.);
  // TH3F* h_histCmTarget = new TH3F("h_histTickness","h_histThickness",NbinX,minHist,1.,NbinX,minHist,1.,NbinY,minHist,1.);
  // TH3F* h_histNameTarget = new TH3F("h_histNameTarget","h_histNameTarget",NbinX,minHist,1.,NbinX,minHist,1.,NbinY,minHist,1.);

  Int_t NbinX = Nbin;
  Int_t NbinY = Nbin;
  double deltaX = (1.-minHist)/static_cast<double>(NbinX);
  double deltaY = (1.-minHist)/static_cast<double>(NbinY);
  
  std::vector<double> Allalpha(NbinX);
  std::vector<double> Allbeta(NbinY);
  for(size_t i = 0;i < Allalpha.size();++i)
    Allalpha[i] = minHist + (i+0.5)*deltaX; 
  for(size_t i = 0;i < Allalpha.size();++i)
    Allbeta[i] = minHist + (i+0.5)*deltaY; 
  
  Histo histo;
  histo.Dohisto = noHisto;
  if(noHisto==false)
    {
      histo.h_hist = new TH2F("h_hist","h_hist",NbinX,minHist,1.,NbinY,minHist,1.);
      histo.h_histNameBeamPrimary = new TH2F("h_histNameBeamPrimary","h_histNameBeamPrimary",NbinX,minHist,1.,NbinY,minHist,1.);
      histo.h_histNameBeamSecondary = new TH2F("h_histNameBeamSecondary","h_histNameBeamSecondary",NbinX,minHist,1.,NbinY,minHist,1.);
      histo.h_histCrossSection = new TH2F("h_histCrosssection","h_histCrossSection",NbinX,minHist,1.,NbinY,minHist,1.);
      histo.h_histEnergy = new TH2F("h_histEnergy","h_histEnergy",NbinX,minHist,1.,NbinY,minHist,1.);
      histo.h_histCmTarget = new TH2F("h_histThickness","h_histThickness",NbinX,minHist,1.,NbinY,minHist,1.);
      histo.h_histNameTarget = new TH2F("h_histNameTarget","h_histNameTarget",NbinX,minHist,1.,NbinY,minHist,1.);
      histo.h_Internal = new TH2F("h_internal","h_internal",15,0,15,1000,-10,10);
    }
  else
    {
      histo.h_hist = nullptr;
      histo.h_histNameBeamPrimary = nullptr;      
      histo.h_histNameBeamSecondary = nullptr;
      histo.h_histCrossSection = nullptr;
      histo.h_histEnergy = nullptr;
      histo.h_histCmTarget = nullptr;
      histo.h_histNameTarget = nullptr;      
      histo.h_Internal = nullptr;
    }
  double maxMax = -99999.;
  std::set<dataMax,CompdataMax> maxMaxId;
  //std::vector<double> maxMaxIdValue(2);
  
  int DisCount = 0;
  //for(int i=1;i<=histo.h_hist->GetNbinsX();++i)
  for(size_t i = 0; i< Allalpha.size();++i)
    {
      
      //for(int j=1;j<=histo.h_hist->GetNbinsY();++j)
      for(size_t j = 0; j<Allbeta.size(); ++j)
	{
	  ++DisCount;
	  if(DisCount==10)
	    {
	      std::cout<<"Loop# ["<<i<<" "<<j<<" ]\n";
	      DisCount=0;
	    }
	  //double alpha = histo.h_hist->GetXaxis()->GetBinCenter(i);
	  //double beta  = histo.h_hist->GetYaxis()->GetBinCenter(j);

	  double alpha = Allalpha[i];
	  double beta  = Allbeta[j];

	  //double gamma  = h_hist->GetZaxis()->GetBinCenter(k);
	  //double gamma = ComputeWeigth(alpha,beta,NtypeComp);

	  if(compute.ParAcceptable(alpha,beta,Ntype) == false)
	    continue;

	  double max_G = -1.e100;
	  double max_CmTarget = -1;
	  std::vector<int> max_nameIds(6,0);
	  double maxIntensity = 0.;
	  double maxContamination = 0.;
	  double maxEnergyM = 0.;
	  double maxEnergyS = 0.;
	  double maxSurvival = 0.;
	  double maxCrossSection = 0;

	  for(auto it_secondaryReaction = hyp_prod.cbegin(), it_secondaryReactionEnd = hyp_prod.cend(); it_secondaryReaction != it_secondaryReactionEnd; ++it_secondaryReaction)
	    {
	      if(withStableBeam==false)
		{
		  auto it_findStable = StableBeam.find(it_secondaryReaction->nameB);
		  if(it_findStable != StableBeam.end())
		    continue;
		}
	      auto it_findTargetS = SecondTargetDensity.find(it_secondaryReaction->nameT);
	      if(it_findTargetS == SecondTargetDensity.end())
		continue;
	      // HypDataProd HypData = *it_secondaryReaction;
	      // //double SecondTargetD = it_findTargetS->second;
	      // HypData.SecondTargetD = it_findTargetS->second;
	      auto HypData = *it_secondaryReaction;
	      //std::cout<<"case:"<<it_secondaryReaction->nameB<<" "<<it_secondaryReaction->nameT<<" CX :["<<it_secondaryReaction->cross_section<<" ]"<<std::endl;
	      for(auto it_FBT = BeamTarget.cbegin(), it_FBT_end = BeamTarget.cend();it_FBT!= it_FBT_end;++it_FBT)
		{
		  auto it_findRejected = Reject.RejectedBeam.find(it_FBT->second.F);
		  if(it_findRejected != Reject.RejectedBeam.end() && (it_FBT->second.BT !="Primary+Primary"))
		    continue;
		  
		  if(HypData.nameB == it_FBT->second.F )
		    {	
		      for(unsigned int idCmTarget = 0; idCmTarget< it_FBT->second.TargetCm.size();++idCmTarget)
			{
			  auto res_Min = minimizer.DoMin(alpha,beta,idCmTarget,it_FBT->second,HypData,histo,Reg);
			  int error_res = std::get<0>(res_Min);
			  if(error_res!=0)
			    continue;
			  double tempG = std::get<1>(res_Min);
			  
			  if(tempG > max_G)
			    {
			      std::tie(std::ignore,max_G,max_CmTarget,max_nameIds,maxIntensity,maxContamination,maxEnergyM,maxEnergyS,maxSurvival,maxCrossSection)=res_Min;

			      //max_G = tempG;
			      //max_CmTarget = t1.TargetCm[idCmTarget];//idCmTarget;
			      //max_nameIds = tempIds;
			      //maxIntensity = t1.ProdFrag[idCmTarget];
			      //maxContamination = t1.ProdPara1[idCmTarget]+t1.ProdPara2[idCmTarget]+t1.ProdPara3[idCmTarget];
			      //maxEnergyM = t1.EnergyMean[idCmTarget];
			      //maxEnergyS = t1.EnergySigma[idCmTarget];
			      //maxSurvival = t1.Survival[idCmTarget];
			      //maxCrossSection = tempICX;//t1.ProdFrag[idCmTarget]*it_secondaryReaction->cross_section;
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
	  
	  if(noHisto==false)
	    {
	      histo.h_hist->Fill(alpha,beta,max_G);
	      histo.h_histNameBeamPrimary->Fill(alpha,beta,nameTemplateB);
	      histo.h_histNameBeamSecondary->Fill(alpha,beta,nameTemplateF);
	      histo.h_histCmTarget->Fill(alpha,beta,max_CmTarget);
	      histo.h_histNameTarget->Fill(alpha,beta,nameTemplateT);
	      histo.h_histEnergy->Fill(alpha,beta,TMath::FloorNint(maxEnergyM*10)*1e-4);
	      histo.h_histCrossSection->Fill(alpha,beta,TMath::FloorNint(maxCrossSection*100)*0.01);
	    }
	  
	}
    }

  if(noHisto==false)
    {
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
      histo.h_Internal->Draw("colz");
      c1->Draw();
      
      TCanvas* c2 = new TCanvas("Cost","Cost",500,500);
      //c2->Divide(2,2);
      // c2->cd(1);
      // h_hist->Project3D("xy")->Draw("colz");
      // c2->cd(2);
      // h_hist->Project3D("xz")->Draw("colz");
      // c2->cd(3);
      // h_hist->Project3D("yz")->Draw("colz");
      //c2->cd(4);
      histo.h_hist->Draw("colz");
      for(unsigned int maxIndex = 0; maxIndex < maxMarkers.size(); ++maxIndex)
	maxMarkers[maxIndex]->Draw("same");
      c2->Draw();
  

      TCanvas* c3 = new TCanvas("NameBeam","NameBeam",500,500);
      //c3->Divide(2,2);
      // c3->cd(1);
      // h_histNameBeam->Project3D("xy")->Draw("colz");  
      // c3->cd(2);
      // h_histNameBeam->Project3D("xz")->Draw("colz");
      // c3->cd(3);
      // h_histNameBeam->Project3D("yz")->Draw("colz");
      //c3->cd(4);
      histo.h_histNameBeamPrimary->Draw(OptionDraw.c_str());
      for(unsigned int maxIndex = 0; maxIndex < maxMarkers.size(); ++maxIndex)
	maxMarkers[maxIndex]->Draw("same");

      c3->Draw();

      TCanvas* c3_1 = new TCanvas("NameTarget","NameTarget",500,500);
      //c3_1->Divide(2,2);
      // c3_1->cd(1);
      // h_histNameTarget->Project3D("xy")->Draw("colz");
      // c3_1->cd(2);
      // h_histNameTarget->Project3D("xz")->Draw("colz");
      // c3_1->cd(3);
      // h_histNameTarget->Project3D("yz")->Draw("colz");
      //c3_1->cd(4);
      histo.h_histNameTarget->Draw(OptionDraw.c_str());
      for(unsigned int maxIndex = 0; maxIndex < maxMarkers.size(); ++maxIndex)
	maxMarkers[maxIndex]->Draw("same");

      c3_1->Draw();

      TCanvas* c3_2 = new TCanvas("NameBeam2","NameBeam2",500,500);
      //c3->Divide(2,2);
      // c3->cd(1);
      // h_histNameBeam->Project3D("xy")->Draw("colz");  
      // c3->cd(2);
      // h_histNameBeam->Project3D("xz")->Draw("colz");
      // c3->cd(3);
      // h_histNameBeam->Project3D("yz")->Draw("colz");
      //c3->cd(4);
      histo.h_histNameBeamSecondary->Draw(OptionDraw.c_str());
      for(unsigned int maxIndex = 0; maxIndex < maxMarkers.size(); ++maxIndex)
	maxMarkers[maxIndex]->Draw("same");

      c3_2->Draw();

      TCanvas* c4 = new TCanvas("Thickness","Thickness",500,500);
      //c4->Divide(2,2);
      //c4->cd(1);
      // h_histCmTarget->Project3D("xy")->Draw("colz");
      // c4->cd(2);
      // h_histCmTarget->Project3D("xz")->Draw("colz");
      // c4->cd(3);
      // h_histCmTarget->Project3D("yz")->Draw("colz");
      //c4->cd(4);
      histo.h_histCmTarget->Draw(OptionDraw.c_str());
      for(unsigned int maxIndex = 0; maxIndex < maxMarkers.size(); ++maxIndex)
	maxMarkers[maxIndex]->Draw("same");

      c4->Draw();

      TCanvas* c5 = new TCanvas("CX","CX",500,500);
      //c3->Divide(2,2);
      // c3->cd(1);
      // h_histNameBeam->Project3D("xy")->Draw("colz");  
      // c3->cd(2);
      // h_histNameBeam->Project3D("xz")->Draw("colz");
      // c3->cd(3);
      // h_histNameBeam->Project3D("yz")->Draw("colz");
      //c3->cd(4);
      histo.h_histCrossSection->Draw(OptionDraw.c_str());
      for(unsigned int maxIndex = 0; maxIndex < maxMarkers.size(); ++maxIndex)
	maxMarkers[maxIndex]->Draw("same");

      c5->Draw();
      
      TCanvas* c6 = new TCanvas("Energy","Energy",500,500);
      //c3->Divide(2,2);
      // c3->cd(1);
      // h_histNameBeam->Project3D("xy")->Draw("colz");  
      // c3->cd(2);
      // h_histNameBeam->Project3D("xz")->Draw("colz");
      // c3->cd(3);
      // h_histNameBeam->Project3D("yz")->Draw("colz");
      //c3->cd(4);
      histo.h_histEnergy->Draw(OptionDraw.c_str());
      for(unsigned int maxIndex = 0; maxIndex < maxMarkers.size(); ++maxIndex)
	maxMarkers[maxIndex]->Draw("same");

      c6->Draw();
    }

  // int maxIndex = 0;
  // for(auto dataId : maxMaxId)
  //   {
  //     std::cout<<" Max ["<<maxIndex<<"] -> ("<<dataId.idI<<", "<<dataId.idJ<<") | alpha="<<dataId.alpha<<" beta="<<dataId.beta<<" "<<std::endl;
  //     std::cout<<" +-> beamP:"<<dataId.beam<<" target:"<<dataId.target<<" IDcm:"<<dataId.cmtarget<<" Int:"<<dataId.intensity<<" | cont:"<<dataId.contamination<<std::endl;
  //     std::cout<<" +-> beamS:"<<dataId.fragment<<" Energy:"<<dataId.energyMean<<" +- "<<dataId.energySigma<<" Survival:"<<dataId.survivalRate<<std::endl;
  //     std::cout<<" +-> CX:"<<dataId.crossSectionHyp<<std::endl;
  //     maxMarkers[maxIndex] = new TMarker(dataId.alpha,dataId.beta,30);
  //     ++maxIndex;
  //   }
    
  return std::make_tuple(maxMaxId.size(),maxMaxId.begin()->beam,maxMaxId.begin()->target,maxMaxId.begin()->cmtarget,maxMaxId.begin()->intensity,maxMaxId.begin()->fragment,maxMaxId.begin()->energyMean,maxMaxId.begin()->crossSectionHyp);
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
  std::set<HypDataProd,CompHypDataProd> hyp_prodAll;

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
	  HypDataProd nameBT;
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

  auto Pred = [&] (const HypDataProd& nameBT) -> bool
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
  
  std::set<HypDataProd,CompHypDataProd> hyp_prod;
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

  Interpole doInterpoling;
  for(std::map<std::string, tuple_graph>::iterator it_FBT = BeamTarget.begin(), it_FBT_end = BeamTarget.end();it_FBT!= it_FBT_end;++it_FBT)
    {
      
      
      if(it_FBT->second.F == "C11") 
	{
	  std::cout<<"*Entry :"<<it_FBT->first<<std::endl;
	  if(/*it_FBT->second.BT_ids[0] == 14 && it_FBT->second.BT_ids[1] == 7 &&*/ it_FBT->second.BT_ids[2] == 9 && it_FBT->second.BT_ids[3] == 4)
	    {
	      //auto temp = doInterpoling.Inter(it_FBT->second);
	      //temp.Print();
	      it_FBT->second.Print();
	    }
	}
    }

  std::cout<<" TEST :"<<std::endl;
  auto it_test = BeamTarget.find("C11Si30Be9");
  auto it_test2 = doInterpoling.Inter(it_test->second);
  it_test2.Print();

  TCanvas* c1= new TCanvas("c1","c1",500,500);
  c1->cd();
  Double_t temp1x [it_test->second.TargetCm.size()];
  Double_t temp1y [it_test->second.TargetCm.size()];
  int tempId = 0;
  for(auto idx : sort_indexes(it_test->second.TargetCm) )
    {
      temp1x[tempId] = it_test->second.TargetCm[idx]; 
      temp1y[tempId] = it_test->second.Survival[idx];
      ++tempId;
    }
  Double_t temp2x [it_test2.TargetCm.size()];
  Double_t temp2y [it_test2.TargetCm.size()];
  tempId = 0;
  for(auto idx : sort_indexes(it_test2.TargetCm) )
    {
      temp2x[tempId] = it_test2.TargetCm[idx]; 
      temp2y[tempId] = it_test2.Survival[idx];
      ++tempId;
    }

  TGraph* graph1 = new TGraph(it_test->second.TargetCm.size(),temp1x,temp1y); 
  TGraph* graph2 = new TGraph(it_test2.TargetCm.size(),temp2x,temp2y); 

  TMultiGraph* GraphM = new TMultiGraph();
  graph1->SetMarkerColor(1);
  graph1->SetMarkerStyle(4);
  GraphM->Add(graph1,"pl");
  graph2->SetMarkerColor(2);
  graph2->SetMarkerStyle(5);
  GraphM->Add(graph2,"p");

  GraphM->Draw("a*");
  
  return 0;

}



int AnalysisAllHyp(std::string namefileEpax, std::string namefileHyp, int style = 0)
{
  
  std::ifstream ifsHyp (namefileHyp.c_str());
  std::set<std::string> hypAll;

  std::string temp_line;
  while(std::getline(ifsHyp,temp_line))
    {
      std::stringstream stream(temp_line);
      std::string hyp;
      stream >> hyp ;
      if(hyp=="#")
	continue;
      hypAll.emplace(hyp);
    }
  std::cout<<"Hypernuclei :";
  for(auto hyp : hypAll)
    std::cout<<" "<<hyp;
  std::cout<<std::endl;

  bool NoHist = true;
  bool WithStable = true;
  
  std::vector<std::string> results_printout;
  const std::string Atom("0123456789");
  for(auto& Hyp : hypAll)
    {
      size_t Ncand; 
      double nameB;
      double nameT;
      double nameF;
      double CmTarget;
      double Intensity;
      double Energy;
      double crossSectionHyp;
      std::tie(Ncand,nameB,nameT,CmTarget,Intensity,nameF,Energy,crossSectionHyp) = AnalysisHyp(namefileEpax,namefileHyp,Hyp,0,"I-C","AB A B",0,0,-1,20,true,WithStable,NoHist);
      if(Ncand>0)
	{
	  int Ap = nameB/100;
	  int Zp = static_cast<int>(nameB)%100;
	  int At = nameT/100;
	  int Zt = static_cast<int>(nameT)%100;
	  int Af = nameF/100;
	  int Zf = static_cast<int>(nameF)%100;

	  std::string SnameF = ElName2[Zf];
	  SnameF+= std::to_string(Af);
	  std::string SnameP = ElName2[Zp];
	  SnameP+=std::to_string(Ap);
	  TString SnameT = ElName2[Zt];
	  SnameT+=std::to_string(At);
	  if(style == 0)
	    {
	      std::cout<<Hyp<<" ["<< Ncand<<"] :"<<SnameP<<"+"<<SnameT<<" with "<<CmTarget<<" cm"<<" and "<<Intensity<<"/s"<<" -> "<<SnameF<<" @ "<<Energy<<" CX="<<crossSectionHyp<<std::endl; 
	    }
	  else if(style == 1)
	    {
	      std::stringstream out;
	      auto it_Namefound = Hyp.find_first_of(Atom);
	      auto it_Namefound_L = Hyp.find_last_of("L");
	      if(it_Namefound != std::string::npos)
		{
		  std::string tempZ = Hyp.substr(0,it_Namefound);
		  std::string tempA = Hyp.substr(it_Namefound,it_Namefound_L);
		  
		  
		  out <<"^{"<<tempA<<"}_/Lambda"<<tempZ<<" & ^{"<<Ap<<"}"<<ElName2[Zp]<<"+ ^{"<<At<<"}"<<ElName2[Zt]<<" & "<<CmTarget<<" & ^{"<<Af<<"}"<<ElName2[Zf]<<" & "<<Energy<<" & "<<Intensity<<" & "<<crossSectionHyp<<" \\ ";  
		  results_printout.push_back(out.str());
		}
	    }
      
	}
    }

  if(style==1)
    {
      for(auto& res : results_printout)
	std::cout<<res<<std::endl;
    }
  
  // TCanvas* c6 = new TCanvas("Energy","Energy",500,500);
  // //c3->Divide(2,2);
  // // c3->cd(1);
  // // h_histNameBeam->Project3D("xy")->Draw("colz");  
  // // c3->cd(2);
  // // h_histNameBeam->Project3D("xz")->Draw("colz");
  // // c3->cd(3);
  // // h_histNameBeam->Project3D("yz")->Draw("colz");
  // //c3->cd(4);
  // histo.h_histEnergy->Draw(OptionDraw.c_str());
  // for(unsigned int maxIndex = 0; maxIndex < maxMarkers.size(); ++maxIndex)
  //   maxMarkers[maxIndex]->Draw("same");



  return 0;
}
