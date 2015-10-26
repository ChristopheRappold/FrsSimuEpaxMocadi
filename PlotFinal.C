#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <list>
#include <string>
#include <cctype>
#include <algorithm>
#include <functional>

#include "TMath.h"
#include "TString.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGraph.h"
#include "TMultiGraph.h"


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


class tuple_graph 
{
public :
  std::string BT;
  std::vector<double> TargetCm;
  std::vector<double> TargetGramCm;
  
  std::vector<double> ProdFrag;
  std::vector<double> ProdPara1;
  std::vector<double> ProdPara2;
  std::vector<double> ProdPara3;

  std::vector< std::vector<double> > ProdStage;
  
  tuple_graph()
  {}
  tuple_graph(const std::string& Name,double TCm, double TgCM, double prodF, double prodP1, double prodP2, double prodP3,const std::vector<double>& stage):
    BT(Name),TargetCm(1,TCm),TargetGramCm(1,TgCM),ProdFrag(1,prodF),ProdPara1(1,prodP1),ProdPara2(1,prodP2),ProdPara3(1,prodP3),ProdStage(1,stage)
  {
    
  }
  ~tuple_graph() {}
  tuple_graph(const tuple_graph& t):BT(t.BT),TargetCm(t.TargetCm),TargetGramCm(t.TargetGramCm),ProdFrag(t.ProdFrag),ProdPara1(t.ProdPara1),ProdPara2(t.ProdPara1),ProdPara3(t.ProdPara3),ProdStage(t.ProdStage)
  {}
  void AddValue(double TCm, double TgCM, double prodF, double prodP1, double prodP2, double prodP3,const std::vector<double>& stage)
  {
    TargetCm.push_back(TCm);
    TargetGramCm.push_back(TgCM);
    ProdFrag.push_back(prodF);
    ProdPara1.push_back(prodP1);
    ProdPara2.push_back(prodP2);
    ProdPara3.push_back(prodP3);
    ProdStage.push_back(stage);
  }
};



void PlotFinal(const std::string& in_file)
{

  std::ifstream ifs ( in_file.c_str() );

  std::string temp_line;
  int Af,Zf,Ap,Zp,At,Zt;
  double TargetDensity,xs,ProdRate;
  int Sister1,Sister2,Sister3,ZZ;
  double CmTarget, gCmTarget;
  int size1;
  std::vector<double> ProdPourcentage;
  std::vector<double> ProdRateFRS;
  int size2;
  std::vector<int> FragID;
  std::vector<double> ProdFinal;


  std::map<std::string,tuple_graph > BeamTarget;

  while(std::getline(ifs,temp_line))
    {
      
      std::stringstream stream(temp_line);
      stream >> Af >> Zf >> Ap >>Zp >> At >> Zt >> TargetDensity >> xs >> ProdRate >> Sister1 >> Sister2 >> Sister3 >> ZZ >> CmTarget >> gCmTarget;
      stream >> size1;
      ProdPourcentage.resize(size1/2);
      ProdRateFRS.resize(size1/2);
      for(int i=0;i<size1/2;++i)
	{
	  double temp1 , temp2;
	  stream >> temp1 >> temp2;

	  ProdPourcentage[i] = temp1;
	  ProdRateFRS[i] = temp2;
	}
      stream >> size2;
      FragID.resize(size2/2);
      ProdFinal.resize(size2/2);
      for(int i=0;i<size2/2;++i)
	{
	  double temp1 , temp2;
	  stream >> temp1 >> temp2;
	  FragID[i] = temp1; 
	  ProdFinal[i] = temp2;
	}

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
	  tuple_graph temp(SnameBT,CmTarget,gCmTarget,ProdFinal[0],ProdFinal[1],ProdFinal[2],ProdFinal[3],ProdRateFRS);
	  BeamTarget.insert(std::pair<std::string,tuple_graph>(SnameAll,temp));
	}
      else
	{
	  it_FBT->second.AddValue(CmTarget,gCmTarget,ProdFinal[0],ProdFinal[1],ProdFinal[2],ProdFinal[3],ProdRateFRS);
	}
    }


  

  TGraph* GraphAll[BeamTarget.size()];
  TGraph* GraphAllRatio[BeamTarget.size()];

  TMultiGraph* GraphM = new TMultiGraph();
  TMultiGraph* GraphRatioM = new TMultiGraph();

  int i = 0;
  for(std::map<std::string, tuple_graph>::const_iterator it_FBT = BeamTarget.begin(), it_FBT_end = BeamTarget.end();it_FBT!= it_FBT_end;++it_FBT)
    {
      //cout<<"Status :"<<it_FBT->first<<" | "<<it_FBT->second.BT<<" "<<it_FBT->second.TargetCm.size()<<endl;

      GraphAll[i] = new TGraph(it_FBT->second.TargetCm.size());
      GraphAllRatio[i] = new TGraph(it_FBT->second.TargetCm.size());

      GraphAll[i]->SetNameTitle(it_FBT->second.BT.c_str(),it_FBT->second.BT.c_str());
      GraphAllRatio[i]->SetNameTitle(it_FBT->second.BT.c_str(),it_FBT->second.BT.c_str());

      for(unsigned int k = 0;k<it_FBT->second.TargetCm.size();++k)
	{
	  GraphAll[i]->SetPoint(k,it_FBT->second.TargetGramCm[k]*1e-3,1e10*it_FBT->second.ProdFrag[k]);
	  double temp_ratio = TMath::Abs(it_FBT->second.ProdPara1[k])>1e-6 ? it_FBT->second.ProdFrag[k]/it_FBT->second.ProdPara1[k] : it_FBT->second.ProdFrag[k]*10e6;
	  GraphAllRatio[i]->SetPoint(k,it_FBT->second.TargetGramCm[k]*1e-3,temp_ratio);
	}
      GraphAll[i]->SetLineColor(i+1);
      GraphAllRatio[i]->SetLineColor(i+1);

      ++i;
      

      //GraphM->Add(GraphAll[i]);
      //GraphRatioM->Add(GraphAllRatio[i]);
    }
  
  for(unsigned int k = 0;k<BeamTarget.size();++k)
    {
      if(0==GraphAll[k])
	cout<<"E> GrapAll ["<<k<<"] NULL !"<<endl;
      else
	GraphM->Add(GraphAll[k]);
      
      if(0==GraphAllRatio[k])
	cout<<"E> GrapAllRatio ["<<k<<"] NULL !"<<endl;
      else
	GraphRatioM->Add(GraphAllRatio[k]);
    }

  //cout<<GraphAll[0]->GetName()<< " "<<GraphAllRatio[0]->GetName()<<endl;
  //cout<<GraphAll[10]->GetName()<< " "<<GraphAllRatio[10]->GetName()<<endl;

  TCanvas* c1 = new TCanvas("c1","c1",500,500);
  c1->cd();
  //GraphAll[0]->Draw("A*L");
  GraphM->Draw("a*l");
  TCanvas* c2 = new TCanvas("c2","c2",500,500);
  c2->cd();
  //GraphAllRatio[0]->Draw("A*L");
  GraphRatioM->Draw("a*l");

}
