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
#include "TRandom3.h"

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



void computeMeanSigma(const std::vector<double>& vec, double* res)
{
  if(vec.size()==1)
    {
      res[0] = vec[0];
      res[1] = 0.000001;
      return;
    }
  else
    {
      double temp_mean = 0.; 
      for(unsigned int i=0;i<vec.size();++i) 
	{ 
	  temp_mean += vec[i];
	}
      temp_mean /= static_cast<double>(vec.size());

      double temp_var = 0;
      for(unsigned int i=0; i<vec.size(); ++i)
	temp_var += (vec[i]-temp_mean)*(vec[i]-temp_mean);
      
      temp_var /= static_cast<double>(vec.size()-1);
      res[0] = temp_mean;
      res[1] = TMath::Sqrt(temp_var);
      return;
    }
}

double ComputeWeigth(double alpha, double beta, int Ntype)
{
  switch(Ntype)
    {
    case 0 :
      return TMath::Sqrt(1. - alpha*alpha - beta*beta) ;
      break;
    case 1 :
      return -TMath::Sqrt(1. - alpha*alpha - beta*beta) ;
      break;
    case 2 :
      return 1. - TMath::Abs(alpha) - TMath::Abs(beta) ;
      break;
    case 3 :
      return TMath::Abs(alpha) + TMath::Abs(beta) - 1.;
    default :
      return 1 - alpha - beta;
      break;
    }
  return gRandom->Uniform(0,1);

  

}

int setPermutation(const std::string& type,const std::vector<double>& in, std::vector<double>& out, int NtypeComp = -1)
{
  if(in.size()!=2 || out.size()!=3)
    {
      std::cout<<"E> Permutation vector in out not valid "<<in.size()<<" "<<out.size()<<std::endl;
      return -1;
    }
  if(type=="A B AB")
    {
      out[0] = in[0];
      out[1] = in[1];
      out[2] = ComputeWeigth(in[0],in[1],NtypeComp);//1-in[0]-in[1];
    }
  else if(type=="B A AB")
    {
      out[0] = in[1];
      out[1] = in[0];
      out[2] = ComputeWeigth(in[0],in[1],NtypeComp);//1-in[0]-in[1];
    }
  else if(type=="A AB B")
    {
      out[0] = in[0];
      out[1] = ComputeWeigth(in[0],in[1],NtypeComp);//1-in[0]-in[1];
      out[2] = in[1];
    }
  else if(type=="B AB A")
    {
      out[0] = in[1];
      out[1] = ComputeWeigth(in[0],in[1],NtypeComp);//1-in[0]-in[1];
      out[2] = in[0];
    }
  else if(type=="AB A B")
    {
      out[0] = ComputeWeigth(in[0],in[1],NtypeComp);//1-in[0]-in[1];
      out[1] = in[0];
      out[2] = in[1];
    }
  else if(type=="AB B A")
    {
      out[0] = ComputeWeigth(in[0],in[1],NtypeComp);//1-in[0]-in[1];
      out[1] = in[1];
      out[2] = in[0];
    }
  else
    {
      std::cout<<"E> Permutation not valid :"<<type<<std::endl;
      return -1;
    }

  return 0;
}

bool ParAcceptable(double alpha, double beta, int Ntype)
{
  switch(Ntype)
    {
    case 0 :
      return (alpha*alpha + beta*beta <= 1.);
      break;
    case 1 :
      return (TMath::Abs(alpha)+TMath::Abs(beta) <= 1.);
      break;
    case 2 :
      return (TMath::Max(TMath::Abs(alpha),TMath::Abs(beta)) <= 1.);
      break;
    default :
      return true;
      break;
    }
  return true;
}


void AnalysisFrag(std::string namefile, int type = 0, std::string ICtype = "I-C", std::string Mtype = "A B AB", int Ntype = -1, int NtypeComp = 1, double minHist = 0)
{

  std::ifstream ifs ( namefile.c_str() );

  std::string temp_line;
  int Af,Zf,Ap,Zp,At,Zt;
  double TargetDensity,xs,ProdRate;
  //int Sister1,Sister2,Sister3,ZZ;
  double CmTarget, gCmTarget;
  int size1;
  int id1;
  double prod1;
  int id2;
  double prod2;
  int id3;
  double prod3;
  int id4;
  double prod4;

  std::map<std::string,tuple_graph > BeamTarget;

  std::vector<double> intensity;
  std::vector<double> contamination;
  std::vector<double> targetLength;

  std::vector<double> IntCon;

  std::getline(ifs,temp_line);

  while(std::getline(ifs,temp_line))
    {
      std::stringstream stream(temp_line);
      stream >> Af >> Zf >> Ap >> Zp >> At >> Zt >> TargetDensity >> xs >> ProdRate >> CmTarget >> gCmTarget >> size1 
	       >> id1 >> prod1 >> id2 >> prod2 >> id3 >> prod3 >> id4 >> prod4;

      intensity.push_back(1e10*prod1);
      contamination.push_back(1e10*prod2+1e10*prod3+1e10*prod4);
      targetLength.push_back(gCmTarget);
      IntCon.push_back(1e10*prod1-(1e10*prod2+1e10*prod3+1e10*prod4));
      std::vector<int> idsBT ;
      idsBT.push_back(Ap);
      idsBT.push_back(Zp);
      idsBT.push_back(At);
      idsBT.push_back(Zt);
      TString nameF("");
      nameF+=Af;
      nameF+= ElName2[Zf];
      TString nameP("");
      nameP+=Ap;
      nameP+= ElName2[Zp];
      TString nameT(""); 
      nameT+=At;
      nameT+= ElName2[Zt];

      std::string SnameF(nameF.Data());
      std::string SnameP(nameP.Data());
      std::string SnameT(nameT.Data());
      std::string SnameAll(SnameF);
      SnameAll+=SnameP;
      SnameAll+=SnameT;
      std::string SnameBT(SnameP);
      SnameBT+="+";
      SnameBT+=SnameT;

      std::map<std::string,tuple_graph>::iterator it_FBT = BeamTarget.find(SnameAll);
      if(it_FBT==BeamTarget.end())
	{
	  tuple_graph temp(SnameBT,idsBT,CmTarget,gCmTarget,1e10*prod1,1e10*prod2,1e10*prod3,1e10*prod4);
	  BeamTarget.insert(std::pair<std::string,tuple_graph>(SnameAll,temp));
	}
      else
	{
	  it_FBT->second.AddValue(CmTarget,gCmTarget,1e10*prod1,1e10*prod2,1e10*prod3,1e10*prod4);
	}

    }

  // for(std::map<std::string, tuple_graph>::const_iterator it_FBT = BeamTarget.begin(), it_FBT_end = BeamTarget.end();it_FBT!= it_FBT_end;++it_FBT)
  //   {
  //     std::cout<<"*Entry :"<<it_FBT->first<<std::endl;
  //     it_FBT->second.Print();
  //   }

  double Tminmax[2] = {*std::min_element(targetLength.begin(),targetLength.end()),*std::max_element(targetLength.begin(),targetLength.end())};
  double Iminmax[2] = {*std::min_element(intensity.begin(),intensity.end()),*std::max_element(intensity.begin(),intensity.end())};
  double Cminmax[2] = {*std::min_element(contamination.begin(),contamination.end()),*std::max_element(contamination.begin(),contamination.end())};
  double ICminmax[2] = {*std::min_element(IntCon.begin(),IntCon.end()),*std::max_element(IntCon.begin(),IntCon.end())};

  double Tmeansigma [2];
  double Imeansigma [2];
  double Cmeansigma [2];
  double ICmeansigma [2];
  
  computeMeanSigma(intensity,Imeansigma);
  computeMeanSigma(targetLength,Tmeansigma);
  computeMeanSigma(contamination,Cmeansigma);
  computeMeanSigma(IntCon,ICmeansigma);

  std::cout<<" mix max "<<Tminmax[0]<<" "<<Tminmax[1]<<" "<<Iminmax[0]<<" "<<Iminmax[1]<<" "<<ICminmax[0]<<" "<<ICminmax[1]<<std::endl;

  std::cout<<" type : "<<type<<" "<<ICtype<<std::endl;
  
  // TH3F* h_hist = new TH3F("h_hist","h_hist",NbinX,minHist,1.,NbinX,minHist,1.,NbinY,minHist,1.);
  // TH3F* h_histNameBeam = new TH3F("h_histNameBeam","h_histNameBeam",NbinX,minHist,1.,NbinX,minHist,1.,NbinY,minHist,1.);
  // TH3F* h_histCmTarget = new TH3F("h_histTickness","h_histThickness",NbinX,minHist,1.,NbinX,minHist,1.,NbinY,minHist,1.);
  // TH3F* h_histNameTarget = new TH3F("h_histNameTarget","h_histNameTarget",NbinX,minHist,1.,NbinX,minHist,1.,NbinY,minHist,1.);

  Int_t NbinX = 200;
  Int_t NbinY = 200;
  TH2F* h_hist = new TH2F("h_hist","h_hist",NbinX,minHist,1.,NbinY,minHist,1.);
  TH2F* h_histNameBeam = new TH2F("h_histNameBeam","h_histNameBeam",NbinX,minHist,1.,NbinY,minHist,1.);
  TH2F* h_histCmTarget = new TH2F("h_histThickness","h_histThickness",NbinX,minHist,1.,NbinY,minHist,1.);
  TH2F* h_histNameTarget = new TH2F("h_histNameTarget","h_histNameTarget",NbinX,minHist,1.,NbinY,minHist,1.);

  TH2F* h_Internal = new TH2F("h_internal","h_internal",10,0,10,100,-4,4);

  double maxMax = -99999.;
  std::vector<dataMax> maxMaxId;
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

	  if(ParAcceptable(alpha,beta,Ntype) == false)
	    continue;

	  double max_G = -1.e100;
	  int max_CmTraget = -1;
	  std::string max_name;
	  std::vector<int> max_nameIds(4,0);
	  double maxIntensity=0.;
	  double maxContamination=0.;
	  for(std::map<std::string, tuple_graph>::const_iterator it_FBT = BeamTarget.begin(), it_FBT_end = BeamTarget.end();it_FBT!= it_FBT_end;++it_FBT)
	    {
	      std::string tempName (it_FBT->second.BT);
	      std::vector<int> tempIds(it_FBT->second.BT_ids);
	      
	      for(unsigned int idCmTarget = 0; idCmTarget< it_FBT->second.TargetCm.size();++idCmTarget)
		{
		  double tempT = it_FBT->second.TargetGramCm[idCmTarget];
		  double tempI = it_FBT->second.ProdFrag[idCmTarget];
		  double tempC = it_FBT->second.ProdPara1[idCmTarget]+it_FBT->second.ProdPara2[idCmTarget]+it_FBT->second.ProdPara3[idCmTarget];
		  if(ICtype == "I-C")
		    tempC = tempI - tempC;

		  h_Internal->Fill(6.,(tempT-Tminmax[0])/(Tminmax[1]-Tminmax[0]));
		  h_Internal->Fill(7.,(tempI-Iminmax[0])/(Iminmax[1]-Iminmax[0]));
		  if(ICtype == "I-C")
		    h_Internal->Fill(8.,(tempC-ICminmax[0])/(ICminmax[1]-ICminmax[0]));
		  else
		    h_Internal->Fill(8.,(tempC-Cminmax[0])/(Cminmax[1]-Cminmax[0]));

		  if(type==0)
		    {
		      tempT -= Tminmax[0];
		      tempT /= Tminmax[1]-Tminmax[0];
		      
		      tempI -= Iminmax[0];
		      tempI /= Iminmax[1]-Iminmax[0];
		      
		      if(ICtype == "I-C")
			{
			  tempC -= ICminmax[0];
			  tempC /= ICminmax[1]-ICminmax[0];
			}
		      else
			{
			  tempC -= Cminmax[0];
			  tempC /= Cminmax[1]-Cminmax[0];
			}
		    }
		  else
		    {
		      tempT -= Tmeansigma[0];
		      tempT /= Tmeansigma[1];

		      tempI -= Imeansigma[0];
		      tempI /= Imeansigma[1];

		      if(ICtype=="I-C")
			{
			  tempC -= ICmeansigma[0];
			  tempC /= ICmeansigma[1];
			}
		      else
			{
			  tempC -= Cmeansigma[0];
			  tempC /= Cmeansigma[1];
			}
		    }
		  // double norm = TMath::Abs(alpha + beta - gamma);
		  // double tempG = alpha*tempI + beta*tempT - gamma*tempC ;
		  // tempG /= norm;
		
		  //double tempG = alpha*tempI - beta*tempT + (1-alpha+beta)*tempC ;
		  std::vector<double> inPar(2);
		  inPar[0] = alpha;
		  inPar[1] = beta;

		  std::vector<double> outPar(3);
		  int res = setPermutation(Mtype,inPar,outPar,NtypeComp);
		  if(res==-1)
		    return;

		  double tempG =  - outPar[0]*tempT + outPar[1]*tempC + outPar[2]*tempI   ;
		  
		  h_Internal->Fill(0.,tempT);
		  h_Internal->Fill(1.,tempI);
		  h_Internal->Fill(2.,tempC);
		  h_Internal->Fill(3.,tempG);
		  
		  if(tempG > max_G)
		    {
		      max_G = tempG;
		      max_name = tempName;
		      max_CmTraget = idCmTarget;
		      max_nameIds = tempIds;
		      maxIntensity = it_FBT->second.ProdFrag[idCmTarget];
		      maxContamination = it_FBT->second.ProdPara1[idCmTarget]+it_FBT->second.ProdPara2[idCmTarget]+it_FBT->second.ProdPara3[idCmTarget];
		    }
		}
	    }
	  
	  
	  double nameTemplateB = max_nameIds[0]*100 + max_nameIds[1];
	  double nameTemplateT = max_nameIds[2]*100 + max_nameIds[3];

	  // h_hist->Fill(alpha,beta,gamma,max_G);
	  // h_histNameBeam->Fill(alpha,beta,gamma,nameTemplateT);
	  // h_histCmTarget->Fill(alpha,beta,gamma,max_CmTraget);
	  // h_histNameTarget->Fill(alpha,beta,gamma,nameTemplateB);
	  if(TMath::Abs(max_G-maxMax)<=1e-4)
	    {
	      
	      dataMax tempMax = {i,j,alpha,beta, nameTemplateB, nameTemplateT, max_CmTraget, maxIntensity, maxContamination};
	      maxMaxId.push_back(tempMax);
	      
	      //std::cout<<" max~ #"<<i<<" "<<j<<" "<<max_G<<" | "<<maxMax <<" size:"<<maxMaxId.size()<<std::endl;
	    }
	  else if(max_G>maxMax)
	    {
	      dataMax tempMax = {i,j, alpha, beta, nameTemplateB, nameTemplateT, max_CmTraget, maxIntensity, maxContamination};
	      //maxMaxId[0] = i;
	      //maxMaxId[1] = j;

	      //maxMaxIdValue[0] = alpha;
	      //maxMaxIdValue[1] = beta;
	      if(maxMaxId.size()==0)
		maxMaxId.push_back(tempMax);
	      else
		{
		  maxMaxId.clear();
		  maxMaxId.push_back(tempMax);
		}
	      //std::cout<<" Max #"<<i<<" "<<j<<" "<<max_G<<" | "<<maxMax<<" size:"<<maxMaxId.size()<<std::endl;
	      
	      maxMax = max_G;
	    }


	  h_hist->Fill(alpha,beta,max_G);
	  h_histNameBeam->Fill(alpha,beta,nameTemplateB);
	  h_histCmTarget->Fill(alpha,beta,max_CmTraget+1);
	  h_histNameTarget->Fill(alpha,beta,nameTemplateT);
	  
	}
    }

  //std::cout<<" Max of the maxG : ["<<maxMaxId[0]<<", "<<maxMaxId[1]<<"] | alpha="<<maxMaxIdValue[0]<<" beta="<<maxMaxIdValue[1]<<std::endl;
  std::cout<<" #max :"<<maxMaxId.size()<<std::endl;
  std::vector<TMarker*> maxMarkers (maxMaxId.size(),0);

  for(unsigned int maxIndex = 0; maxIndex < maxMaxId.size(); ++maxIndex)
    {
      std::cout<<" Max ["<<maxIndex<<"] -> ("<<maxMaxId[maxIndex].idI<<", "<<maxMaxId[maxIndex].idJ<<") | alpha="<<maxMaxId[maxIndex].alpha<<" beta="<<maxMaxId[maxIndex].beta<<" "<<std::endl;
      std::cout<<" +-> beam:"<<maxMaxId[maxIndex].beam<<" target:"<<maxMaxId[maxIndex].target<<" IDcm:"<<maxMaxId[maxIndex].cmtarget<<" Int:"<<maxMaxId[maxIndex].intensity<<" | cont:"<<maxMaxId[maxIndex].contamination<<std::endl;
      maxMarkers[maxIndex] = new TMarker(maxMaxId[maxIndex].alpha,maxMaxId[maxIndex].beta,30);
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
  h_histNameBeam->Draw(OptionDraw.c_str());
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

}


int AnalysisHyp(std::string namefileEpax, std::string namefileHyp, std::string Hyp, int type = 0, std::string ICtype = "I-C", std::string Mtype = "A B AB", int Ntype = -1, int NtypeComp = 1, double minHist = 0)
{
  double PrimaryBeamInt = 1e10;

  std::ifstream ifsHyp (namefileHyp.c_str());
  std::string temp_line;
  std::set<nameBTtuple,CompBTTuple> hyp_prod;

  std::getline(ifsHyp,temp_line);
    
  while(std::getline(ifsHyp,temp_line))
    {
      std::stringstream stream(temp_line);
      std::string hyp;
      int nb;
      stream >> hyp >> nb;
      if(hyp != Hyp)
	continue;

      nb /= 2;
      for(int id = 0; id<nb ; ++id)
	{
	  std::string BT;
	  double CX1;

	  stream >> BT >> CX1;
	  
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

	  hyp_prod.insert(nameBT);
	  
	}
    }

  
  
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
  std::vector<double> EnergyM;
  std::vector<double> Survival;
  std::vector<double> Trans;
  std::vector<double> IntCon;

  std::getline(ifs,temp_line);

  while(std::getline(ifs,temp_line))
    {
      std::stringstream stream(temp_line);
      stream >> Af >> Zf >> Ap >> Zp >> At >> Zt >> TargetDensity >> xs >> ProdRate >> CmTarget >> Sister1 >> Sister2 >> Sister3 >> ZZ
	     >> CmTarget >> SurvivalRate >> MeanEnergy >> SigmaEnergy >> Transmission >> prod1 >> prod2 >> prod3 >> prod4 ; 

      
      intensity.push_back(PrimaryBeamInt*prod1);
      contamination.push_back(PrimaryBeamInt*prod2+PrimaryBeamInt*prod3+PrimaryBeamInt*prod4);
      targetLength.push_back(gCmTarget);
      IntCon.push_back(PrimaryBeamInt*prod1-(PrimaryBeamInt*prod2+PrimaryBeamInt*prod3+PrimaryBeamInt*prod4));
      
      std::vector<int> idsBT ;
      idsBT.push_back(Ap);
      idsBT.push_back(Zp);
      idsBT.push_back(At);
      idsBT.push_back(Zt);
      TString nameF("");
      nameF+=Af;
      nameF+= ElName2[Zf];
      TString nameP("");
      nameP+=Ap;
      nameP+= ElName2[Zp];
      TString nameT(""); 
      nameT+=At;
      nameT+= ElName2[Zt];

      std::string SnameF(nameF.Data());
      std::string SnameP(nameP.Data());
      std::string SnameT(nameT.Data());
      std::string SnameAll(SnameF);
      SnameAll+=SnameP;
      SnameAll+=SnameT;
      std::string SnameBT(SnameP);
      SnameBT+="+";
      SnameBT+=SnameT;

      std::map<std::string,tuple_graph>::iterator it_FBT = BeamTarget.find(SnameAll);
      if(it_FBT==BeamTarget.end())
	{
	  tuple_graph temp(SnameBT,idsBT,CmTarget,gCmTarget,PrimaryBeamInt*prod1,PrimaryBeamInt*prod2,PrimaryBeamInt*prod3,PrimaryBeamInt*prod4);
	  BeamTarget.insert(std::pair<std::string,tuple_graph>(SnameAll,temp));
	}
      else
	{
	  it_FBT->second.AddValue(CmTarget,gCmTarget,PrimaryBeamInt*prod1,PrimaryBeamInt*prod2,PrimaryBeamInt*prod3,PrimaryBeamInt*prod4);
	}

    }

  // for(std::map<std::string, tuple_graph>::const_iterator it_FBT = BeamTarget.begin(), it_FBT_end = BeamTarget.end();it_FBT!= it_FBT_end;++it_FBT)
  //   {
  //     std::cout<<"*Entry :"<<it_FBT->first<<std::endl;
  //     it_FBT->second.Print();
  //   }

  double Tminmax[2] = {*std::min_element(targetLength.begin(),targetLength.end()),*std::max_element(targetLength.begin(),targetLength.end())};
  double Iminmax[2] = {*std::min_element(intensity.begin(),intensity.end()),*std::max_element(intensity.begin(),intensity.end())};
  double Cminmax[2] = {*std::min_element(contamination.begin(),contamination.end()),*std::max_element(contamination.begin(),contamination.end())};
  double ICminmax[2] = {*std::min_element(IntCon.begin(),IntCon.end()),*std::max_element(IntCon.begin(),IntCon.end())};

  double Tmeansigma [2];
  double Imeansigma [2];
  double Cmeansigma [2];
  double ICmeansigma [2];
  
  computeMeanSigma(intensity,Imeansigma);
  computeMeanSigma(targetLength,Tmeansigma);
  computeMeanSigma(contamination,Cmeansigma);
  computeMeanSigma(IntCon,ICmeansigma);

  std::cout<<" mix max "<<Tminmax[0]<<" "<<Tminmax[1]<<" "<<Iminmax[0]<<" "<<Iminmax[1]<<" "<<ICminmax[0]<<" "<<ICminmax[1]<<std::endl;

  std::cout<<" type : "<<type<<" "<<ICtype<<std::endl;
  
  // TH3F* h_hist = new TH3F("h_hist","h_hist",NbinX,minHist,1.,NbinX,minHist,1.,NbinY,minHist,1.);
  // TH3F* h_histNameBeam = new TH3F("h_histNameBeam","h_histNameBeam",NbinX,minHist,1.,NbinX,minHist,1.,NbinY,minHist,1.);
  // TH3F* h_histCmTarget = new TH3F("h_histTickness","h_histThickness",NbinX,minHist,1.,NbinX,minHist,1.,NbinY,minHist,1.);
  // TH3F* h_histNameTarget = new TH3F("h_histNameTarget","h_histNameTarget",NbinX,minHist,1.,NbinX,minHist,1.,NbinY,minHist,1.);

  Int_t NbinX = 200;
  Int_t NbinY = 200;
  TH2F* h_hist = new TH2F("h_hist","h_hist",NbinX,minHist,1.,NbinY,minHist,1.);
  TH2F* h_histNameBeam = new TH2F("h_histNameBeam","h_histNameBeam",NbinX,minHist,1.,NbinY,minHist,1.);
  TH2F* h_histCmTarget = new TH2F("h_histThickness","h_histThickness",NbinX,minHist,1.,NbinY,minHist,1.);
  TH2F* h_histNameTarget = new TH2F("h_histNameTarget","h_histNameTarget",NbinX,minHist,1.,NbinY,minHist,1.);

  TH2F* h_Internal = new TH2F("h_internal","h_internal",10,0,10,100,-4,4);

  double maxMax = -99999.;
  std::vector<dataMax> maxMaxId;
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

	  if(ParAcceptable(alpha,beta,Ntype) == false)
	    continue;

	  double max_G = -1.e100;
	  int max_CmTraget = -1;
	  std::string max_name;
	  std::vector<int> max_nameIds(4,0);
	  double maxIntensity=0.;
	  double maxContamination=0.;
	  for(std::map<std::string, tuple_graph>::const_iterator it_FBT = BeamTarget.begin(), it_FBT_end = BeamTarget.end();it_FBT!= it_FBT_end;++it_FBT)
	    {
	      std::string tempName (it_FBT->second.BT);
	      std::vector<int> tempIds(it_FBT->second.BT_ids);
	      
	      for(unsigned int idCmTarget = 0; idCmTarget< it_FBT->second.TargetCm.size();++idCmTarget)
		{
		  double tempT = it_FBT->second.TargetGramCm[idCmTarget];
		  double tempI = it_FBT->second.ProdFrag[idCmTarget];
		  double tempC = it_FBT->second.ProdPara1[idCmTarget]+it_FBT->second.ProdPara2[idCmTarget]+it_FBT->second.ProdPara3[idCmTarget];
		  if(ICtype == "I-C")
		    tempC = tempI - tempC;

		  h_Internal->Fill(6.,(tempT-Tminmax[0])/(Tminmax[1]-Tminmax[0]));
		  h_Internal->Fill(7.,(tempI-Iminmax[0])/(Iminmax[1]-Iminmax[0]));
		  if(ICtype == "I-C")
		    h_Internal->Fill(8.,(tempC-ICminmax[0])/(ICminmax[1]-ICminmax[0]));
		  else
		    h_Internal->Fill(8.,(tempC-Cminmax[0])/(Cminmax[1]-Cminmax[0]));

		  if(type==0)
		    {
		      tempT -= Tminmax[0];
		      tempT /= Tminmax[1]-Tminmax[0];
		      
		      tempI -= Iminmax[0];
		      tempI /= Iminmax[1]-Iminmax[0];
		      
		      if(ICtype == "I-C")
			{
			  tempC -= ICminmax[0];
			  tempC /= ICminmax[1]-ICminmax[0];
			}
		      else
			{
			  tempC -= Cminmax[0];
			  tempC /= Cminmax[1]-Cminmax[0];
			}
		    }
		  else
		    {
		      tempT -= Tmeansigma[0];
		      tempT /= Tmeansigma[1];

		      tempI -= Imeansigma[0];
		      tempI /= Imeansigma[1];

		      if(ICtype=="I-C")
			{
			  tempC -= ICmeansigma[0];
			  tempC /= ICmeansigma[1];
			}
		      else
			{
			  tempC -= Cmeansigma[0];
			  tempC /= Cmeansigma[1];
			}
		    }
		  // double norm = TMath::Abs(alpha + beta - gamma);
		  // double tempG = alpha*tempI + beta*tempT - gamma*tempC ;
		  // tempG /= norm;
		
		  //double tempG = alpha*tempI - beta*tempT + (1-alpha+beta)*tempC ;
		  std::vector<double> inPar(2);
		  inPar[0] = alpha;
		  inPar[1] = beta;

		  std::vector<double> outPar(3);
		  int res = setPermutation(Mtype,inPar,outPar,NtypeComp);
		  if(res==-1)
		    return -1;

		  double tempG =  - outPar[0]*tempT + outPar[1]*tempC + outPar[2]*tempI   ;
		  
		  h_Internal->Fill(0.,tempT);
		  h_Internal->Fill(1.,tempI);
		  h_Internal->Fill(2.,tempC);
		  h_Internal->Fill(3.,tempG);
		  
		  if(tempG > max_G)
		    {
		      max_G = tempG;
		      max_name = tempName;
		      max_CmTraget = idCmTarget;
		      max_nameIds = tempIds;
		      maxIntensity = it_FBT->second.ProdFrag[idCmTarget];
		      maxContamination = it_FBT->second.ProdPara1[idCmTarget]+it_FBT->second.ProdPara2[idCmTarget]+it_FBT->second.ProdPara3[idCmTarget];
		    }
		}
	    }
	  
	  
	  double nameTemplateB = max_nameIds[0]*100 + max_nameIds[1];
	  double nameTemplateT = max_nameIds[2]*100 + max_nameIds[3];

	  // h_hist->Fill(alpha,beta,gamma,max_G);
	  // h_histNameBeam->Fill(alpha,beta,gamma,nameTemplateT);
	  // h_histCmTarget->Fill(alpha,beta,gamma,max_CmTraget);
	  // h_histNameTarget->Fill(alpha,beta,gamma,nameTemplateB);
	  if(TMath::Abs(max_G-maxMax)<=1e-4)
	    {
	      
	      dataMax tempMax = {i,j,alpha,beta, nameTemplateB, nameTemplateT, max_CmTraget, maxIntensity, maxContamination};
	      maxMaxId.push_back(tempMax);
	      
	      //std::cout<<" max~ #"<<i<<" "<<j<<" "<<max_G<<" | "<<maxMax <<" size:"<<maxMaxId.size()<<std::endl;
	    }
	  else if(max_G>maxMax)
	    {
	      dataMax tempMax = {i,j, alpha, beta, nameTemplateB, nameTemplateT, max_CmTraget, maxIntensity, maxContamination};
	      //maxMaxId[0] = i;
	      //maxMaxId[1] = j;

	      //maxMaxIdValue[0] = alpha;
	      //maxMaxIdValue[1] = beta;
	      if(maxMaxId.size()==0)
		maxMaxId.push_back(tempMax);
	      else
		{
		  maxMaxId.clear();
		  maxMaxId.push_back(tempMax);
		}
	      //std::cout<<" Max #"<<i<<" "<<j<<" "<<max_G<<" | "<<maxMax<<" size:"<<maxMaxId.size()<<std::endl;
	      
	      maxMax = max_G;
	    }


	  h_hist->Fill(alpha,beta,max_G);
	  h_histNameBeam->Fill(alpha,beta,nameTemplateB);
	  h_histCmTarget->Fill(alpha,beta,max_CmTraget+1);
	  h_histNameTarget->Fill(alpha,beta,nameTemplateT);
	  
	}
    }

  //std::cout<<" Max of the maxG : ["<<maxMaxId[0]<<", "<<maxMaxId[1]<<"] | alpha="<<maxMaxIdValue[0]<<" beta="<<maxMaxIdValue[1]<<std::endl;
  std::cout<<" #max :"<<maxMaxId.size()<<std::endl;
  std::vector<TMarker*> maxMarkers (maxMaxId.size(),0);

  for(unsigned int maxIndex = 0; maxIndex < maxMaxId.size(); ++maxIndex)
    {
      std::cout<<" Max ["<<maxIndex<<"] -> ("<<maxMaxId[maxIndex].idI<<", "<<maxMaxId[maxIndex].idJ<<") | alpha="<<maxMaxId[maxIndex].alpha<<" beta="<<maxMaxId[maxIndex].beta<<" "<<std::endl;
      std::cout<<" +-> beam:"<<maxMaxId[maxIndex].beam<<" target:"<<maxMaxId[maxIndex].target<<" IDcm:"<<maxMaxId[maxIndex].cmtarget<<" Int:"<<maxMaxId[maxIndex].intensity<<" | cont:"<<maxMaxId[maxIndex].contamination<<std::endl;
      maxMarkers[maxIndex] = new TMarker(maxMaxId[maxIndex].alpha,maxMaxId[maxIndex].beta,30);
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
  h_histNameBeam->Draw(OptionDraw.c_str());
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

}


