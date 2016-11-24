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
#include <set>
#include <assert.h>

#include "TMath.h"
#include "TString.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TGeoManager.h"
#include "TGeoElement.h"
#include "TStyle.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TRandom3.h"


static const std::string ElName[110] = {
          "H ","He","Li","Be","B ","C ","N ","O ","F ","Ne","Na","Mg",
          "Al","Si","P ","S ","Cl","Ar","K ","Ca","Sc","Ti","V ","Cr",
          "Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
          "Rb","Sr","Y ","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd",
          "In","Sn","Sb","Te","I ","Xe","Cs","Ba","La","Ce","Pr","Nd",
          "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf",
          "Ta","W ","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po",
          "At","Rn","Fr","Ra","Ac","Th","Pa","U ","Np","Pu","Am","Cm",
          "Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs",
          "Mt","Ds" };


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


using namespace std;

void Extract(std::string name_in,int Atest=9,int Ztest=6,bool updated=false)
{

  std::ifstream ifss("");


  std::ifstream ifs ( name_in.c_str() );

  std::string temp_line;
  int Af,Zf,Ap,Zp,At,Zt;
  double xs,DeXs;
  int sister,ZtZp;

  std::getline(ifs,temp_line);
  
  
  //std::map<int,std::map<int,std::map<int,double> > > table_beam;
  //int beam_passed=0;
  //int beam_all=0;
  double max_xs = -1;
  TString name_max[2];
  cout<<ElName[Ztest-1]<<Atest<<endl;
  TH2F* h_F = new TH2F("BeamTarget","BeamTarget",200,0,200,200,0,200);
  while(std::getline(ifs,temp_line))
    {
      
      std::stringstream stream(temp_line);
      stream >> Af >> Zf >> Ap >>Zp >> At >> Zt >> xs;
      if(updated)
	stream >> DeXs >> sister >> ZtZp ;

      if(Ap==Atest && Zp==Ztest)
	continue;
      TString nameF;
      nameF = Zf==0 ? "n" : ElName[Zf-1];
      TString nameP;
      nameP = Zp==0 ? "n" : ElName[Zp-1];
      TString nameT;
      nameT= Zt ==0 ? "n" : ElName[Zt-1];
    
      nameF+=Af;
      nameP+=Ap;
      nameT+=At;
      
      if(Af==Atest && Zf==Ztest)
	{
	  cout<<nameF<<" "<<nameP<<" "<<nameT<<" "<<DeXs<<endl;
	  if(!updated)
	    {
	      if(xs>max_xs)
		{
		  name_max[0]=nameP;
		  name_max[1]=nameT;
		  max_xs = xs;
		  
		}
	      h_F->Fill(nameP,nameT,xs);
	    }
	  else
	    {
	      if(DeXs>max_xs)
		{
		  name_max[0]=nameP;
		  name_max[1]=nameT;
		  max_xs = DeXs;
		  
		}
	      h_F->Fill(nameP,nameT,DeXs);
	    }
	}
    }
      //std::cout<<beam_passed<<" / "<<beam_all<<endl;
  h_F->Draw("colz text");
  h_F->SetXTitle("Beam");
  h_F->SetYTitle("Target");

  cout<<"Max at : Beam "<<name_max[0]<<" Target "<<name_max[1]<<" Xs"<<max_xs<<endl;
}

void ExtractAll(std::string name_in,bool updated=false, double EkI= -1. /*GeV*/)
{
  //gStyle->SetCanvasPreferGL(kTRUE);

  //std::ifstream ifss("");
  TGeoManager *geom = new TGeoManager("geom","radionuclides");
  TGeoElementTable *tableRN = geom->GetElementTable();

  
  std::ifstream ifs ( name_in.c_str() );

  std::string temp_line;
  int Af,Zf,Ap,Zp,At,Zt;
  double xs,DeXs;
  int sister,ZtZp;

  //std::getline(ifs,temp_line);
  
  std::vector<std::vector<double> > table_XS(90,std::vector<double>(90,-1.));
  std::vector<std::vector<double> > table_XSAll(90,std::vector<double>(90,-1.));
  std::vector<std::vector<TString> > table_Name(90,std::vector<TString>(90));
  std::vector<std::vector<TString> > table_NameAll(90,std::vector<TString>(90));

  std::map<std::string,std::vector<std::string> > BeamTarget;
  //std::map<int,std::map<int,std::map<std::string,double> > > table_BT;
  //int beam_passed=0;
  //int beam_all=0;
  //double max_xs = -1;
  TString name_max[2];
  //cout<<ElName[Ztest-1]<<Atest<<endl;
  //TH2F* h_F = new TH2F("BeamTarget","BeamTarget",50,0,50,50,0,50);
  TH2F* h_Fxs = new TH2F("BeamTarget_XS","BeamTarget_XS",90,0,90,90,0,90);
  TH2F* h_Fxs_Frs = new TH2F("BeamTarget_XS_Frs","BeamTarget_XS_Frs",90,0,90,90,0,90);
  TH2F* h_Fxs_SuperFrs = new TH2F("BeamTarget_XS_SuperFrs","BeamTarget_XS_SuperFrs",90,0,90,90,0,90);
  
  double u = 0.931494061; // GeV
  double m_eminus = 5.1099891e-04; // GeV
  //int current_Z=-1,current_A=-1;

  while(std::getline(ifs,temp_line))
    {
      
      std::stringstream stream(temp_line);
      stream >> Af >> Zf >> Ap >>Zp >> At >> Zt >> xs;
      if(updated)
	stream >> DeXs >> sister >> ZtZp ;
      // if(current_Z==-1 && current_A==-1)
      // 	{
      // 	  current_Z=Zf;
      // 	  current_A=Af;
      // 	}
      // if(current_Z!=Zf && current_A!=Af)
      // 	{
      // 	  TString temp_name2(name_max[0]);
      // 	  temp_name2+="+";
      // 	  temp_name2+=name_max[1];
      // 	  std::string temp_nameTP(temp_name2.Data());


      // 	  table_BT[current_Z][current_A][temp_nameTP]=max_xs;
      // 	  max_xs=-1;
      // 	  current_Z=Zf;
      // 	  current_A=Af;
      // 	}

      TString nameF("^{");
      nameF+=Af;
      nameF+="}";
      nameF+= Zf==0 ? "n" : ElName[Zf-1];
      TString nameP("^{");
      nameP+=Ap;
      nameP+="}";
      nameP+= Zp==0 ? "n" : ElName[Zp-1];
      TString nameT("^{"); 
      nameT+=At;
      nameT+="}";
      nameT+= Zt ==0 ? "n" : ElName[Zt-1];

      TString temp_nameBeamTarget(nameP);
      temp_nameBeamTarget+="+";
      temp_nameBeamTarget+=nameT;
	 
      std::string nameBeamTarget(temp_nameBeamTarget.Data());
      std::map<std::string,std::vector<std::string> >::iterator it_BT = BeamTarget.find(nameBeamTarget);
      if(it_BT==BeamTarget.end())
	{
	  std::string SnameF(nameF.Data());
	  std::vector<std::string> temp_vec(1,SnameF);
	  BeamTarget.insert(std::pair<std::string,std::vector<std::string> >(nameBeamTarget,temp_vec));
	}
      else
	{
	  std::string SnameF(nameF.Data());
	  it_BT->second.push_back(SnameF);
	}
	

      if(Zf!= Zp || Af!=Ap)
	{
	  if(!updated)
	    {
	      if(xs>table_XS[Zf][Af])
		{
		  name_max[0]=nameP;
		  name_max[1]=nameT;
		  TString temp_name2(name_max[0]);
		  temp_name2+="+";
		  temp_name2+=name_max[1];
		  //temp_name2+="\n";
		//temp_name2+=nameF;
		  table_XS[Zf][Af] = xs;
		  table_Name[Zf][Af] = temp_name2;
		}
	    }
	  else
	    {
	      if(DeXs>table_XS[Zf][Af])
		{
		  name_max[0]=nameP;
		  name_max[1]=nameT;
		  TString temp_name2(name_max[0]);
		  temp_name2+="+";
		  temp_name2+=name_max[1];
		  //temp_name2+="\n";
		  //temp_name2+=nameF;
		  table_XS[Zf][Af] = DeXs;
		  table_Name[Zf][Af] = temp_name2;
		}
	      
	    }
	}
      if(!updated)
	{
	  if(xs>table_XSAll[Zf][Af] )
	    {
	      name_max[0]=nameP;
	      name_max[1]=nameT;
	      TString temp_name2(name_max[0]);
	      temp_name2+="+";
	      temp_name2+=name_max[1];
	      //temp_name2+="\n";
	      //temp_name2+=nameF;
	      table_XSAll[Zf][Af] = xs;
	      table_NameAll[Zf][Af] = temp_name2;
	    }
	}
      else
	{
	  if(DeXs>table_XSAll[Zf][Af] )
	    {
	      name_max[0]=nameP;
	      name_max[1]=nameT;
	      TString temp_name2(name_max[0]);
	      temp_name2+="+";
	      temp_name2+=name_max[1];
	      //temp_name2+="\n";
	      //temp_name2+=nameF;
	      table_XSAll[Zf][Af] = DeXs;
	      table_NameAll[Zf][Af] = temp_name2;
	    }
	}

    }
      //std::cout<<beam_passed<<" / "<<beam_all<<endl;
  
  TCanvas* c1 = new TCanvas("c1","c1");
  c1->cd(1);
  c1->SetBorderSize(0);
  c1->SetFillColor(0);
  Int_t Xmin = 0;
  Int_t Xmax = 70;
  Int_t Zmin = 0;
  //Int_t Zmax = 50;
  c1->Range(Xmin-1,Zmin-1,Xmax+2,(Xmax-Xmin+3)/sqrt(2)+Zmin-1);
  //h_F->Draw("text");
  //h_F->SetXTitle("N");
  //h_F->SetYTitle("Z");


  //  for(std::map<int,std::map<int,std::map<std::string,double> > >::const_iterator it_Z = table_BT.begin(),it_Z_end=table_BT.end();it_Z!=it_Z_end;++it_Z)
  //for(std::map<int,std::map<std::string,double> >::const_iterator it_A = it_Z->second.begin(),it_A_end=it_Z->second.end();it_A!=it_A_end;++it_A)
  //   for(std::map<std::string,double>::const_iterator it_BT = it_A->second.begin(),it_BT_end=it_A->second.end();it_BT!=it_BT_end;++it_BT)
  for(unsigned int Z = 0;Z<table_XS.size();++Z)
    for(unsigned int A = 0;A<table_XS[Z].size();++A)
      if(table_XS[Z][A]>0)
	{
	  
	  TString nameF("^{");
	  nameF+=A;
	  nameF+="}";
	  nameF+= Z==0 ? "n" : ElName[Z-1];

	  TPaveText  *b;
	  Int_t N = A - Z ;
	  b = new TPaveText(2*N,2*Z,2*N+2,2*Z+2);
	  b->AddText(table_Name[Z][A].Data());
	  b->AddText(nameF.Data());
	  b->SetFillStyle(1001);
	  //b->SetTextFont(42);
	  //b->SetTextSize(.01);
	  b->SetTextAlign(22);
	  b->SetLineStyle(1);
	  b->SetLineColor(1);
	  b->SetBorderSize(1);
	  b->SetTextAngle(45);
	  b->Draw();
	  
	  //h_F->Fill(it_A->first-it_Z->first,it_Z->first,it_BT->first.c_str());
	  h_Fxs->Fill(N,Z);//,table_XS[Z][A]);
	  if(EkI>0.)
	    {
	      auto* TempElement = tableRN->GetElementRN(A,Z);
	      double Dmass = 0.;
	      if(TempElement==nullptr)
		{
		  cout<<"E> no element ! "<<A<<" "<<Z<<" "<<TempElement<<endl;
		  if(A==23 && Z==14)
		    Dmass = 23.073*1e-3;
		  else
		    continue;
		}
	      else
		Dmass = TempElement->MassEx()*1e-3;
	      double Mass = Dmass+A*u - Z*m_eminus;
	      //double Energy = EkI*A + mass;
	      double p = TMath::Sqrt((EkI*A+Mass)*(EkI*A+Mass) - Mass*Mass);
	      double Brho = 3.10715497*p/Z;
	      if(Brho<18.0)
		h_Fxs_Frs->Fill(N,Z);//,table_XS[Z][A]);
	      if(Brho<20.0)
		h_Fxs_SuperFrs->Fill(N,Z);//,table_XS[Z][A]);
	      
	    }
	}

  TCanvas* c1_1 = new TCanvas("c1_1","c1_1");
  c1_1->cd(1);
  c1_1->SetBorderSize(0);
  c1_1->SetFillColor(0);
  //Int_t Zmax = 50;
  c1_1->Range(Xmin-1,Zmin-1,Xmax+2,(Xmax-Xmin+3)/sqrt(2)+Zmin-1);

  for(unsigned int Z = 0;Z<table_XSAll.size();++Z)
    for(unsigned int A = 0;A<table_XSAll[Z].size();++A)
      if(table_XSAll[Z][A]>0)
	{
	  
	  TString nameF("^{");
	  nameF+=A;
	  nameF+="}";
	  nameF+= Z==0 ? "n" : ElName[Z-1];

	  TPaveText  *b;
	  Int_t N = A - Z ;
	  b = new TPaveText(2.5*N,2.5*Z,2.5*N+2.5,2.5*Z+2.5);
	  b->AddText(table_NameAll[Z][A].Data());
	  b->AddText(nameF.Data());
	  b->SetFillStyle(1001);
	  //b->SetTextFont(42);
	  //b->SetTextSize(.01);
	  b->SetTextAlign(22);
	  b->SetLineStyle(1);
	  b->SetLineColor(1);
	  b->SetBorderSize(1);
	  b->SetTextAngle(45);
	  b->Draw();
	  
	  //h_F->Fill(it_A->first-it_Z->first,it_Z->first,it_BT->first.c_str());
	}


  TCanvas* c2 = new TCanvas("c2","c2");
  c2->cd(1);
  h_Fxs->Draw("box");
  if(EkI>0)
    {
      h_Fxs_SuperFrs->SetFillColor(4);
      h_Fxs_SuperFrs->Draw("box1 same");
      h_Fxs_Frs->SetFillColor(2);
      h_Fxs_Frs->Draw("box1 same");
    }
  h_Fxs->SetXTitle("N");
  h_Fxs->SetYTitle("Z");



  ofstream outfile;
  std::string outname(name_in);
  std::string ext("result.tex");
  outname+=ext;
  outfile.open(outname.c_str());
  outfile<<"\\documentclass[a4paper,11pt]{article}"<<endl;
  outfile<<"\\usepackage{breqn}"<<endl;
  outfile<<"\\begin{document}"<<endl;
  
  for(std::map<std::string,std::vector<std::string> >::const_iterator it_BT=BeamTarget.begin(),it_BT_end=BeamTarget.end();it_BT!=it_BT_end;++it_BT)
    {
      outfile<<"\\begin{dmath} "<<it_BT->first<<" ["<<it_BT->second.size()<<"] : ";
      for(unsigned int i=0;i<it_BT->second.size();++i)
	{
	  outfile<<it_BT->second[i]<<"\\ ,\\ ";
      	}
      outfile<<"\\end{dmath}"<<endl;
      outfile<<endl;
    }
  
  outfile<<"\\end{document}"<<endl;
  outfile.close();
}

class TupleElement 
{
public :
  std::string name;
  int A;
  int Z;
  double xs;
  
  TupleElement():name(""),A(0),Z(0),xs(-1)
  {}
  TupleElement(std::string N,int Af, int Zf, double XS):name(N),A(Af),Z(Zf),xs(XS)
  {}
  ~TupleElement() {}
  TupleElement(const TupleElement& t):name(t.name),A(t.A),Z(t.Z),xs(t.xs)
  {}
};

int TestBeam(double Th=0.75)
{
  std::ifstream ifss("ListBeam.dat");
  std::string temp_lineB,nameBT;
  int As,Zs;
  double density;
  map<double,std::string> table;
  while(std::getline(ifss,temp_lineB))
    {
      std::stringstream stream(temp_lineB);
      stream >> nameBT >> As >> Zs >> density;
      double dd = 1e-27*density/As*6.02214129e23;
      if(Zs>20)
	table.insert(std::pair<double, std::string>(dd,nameBT));
      //cout<<nameBT<<" "<<1e-27*density/As*6.02214129e23<<endl;
    }
  auto it_table_last = table.rbegin();
  double Dmax = it_table_last->first*Th;
  
  for( auto it_table = table.begin(), it_table_end = table.end(); it_table!=it_table_end;++it_table)
    {
      if(it_table->first>Dmax)
	{
	  cout<<it_table->second<<" "<<it_table->first<<endl;
	}
    }
  return 1;
}

int UpdateFile(std::string name_in,std::string option="all",const std::vector<std::string>& vecFrag = std::vector<std::string>())
{
  bool not_all = false;
  if(option!="all")
    {
      not_all = true;
    }
  
  std::ifstream ifss("ListBeam.dat");
    

  //std::vector<std::vector<bool> > table_stable(90,std::vector<bool>(90,false));
  std::vector<std::vector<double> > table_stable(90,std::vector<double>(90,-1.));
  std::string temp_lineB,nameBT;
  int As,Zs;
  double density;
  while(std::getline(ifss,temp_lineB))
    {
      std::stringstream stream(temp_lineB);
      stream >> nameBT >> As >> Zs >> density;
      table_stable[As][Zs]=density;
    }  

  std::ifstream ifss2("ListTarget.dat");
    

  //std::vector<std::vector<bool> > table_stable(90,std::vector<bool>(90,false));
  std::vector<std::vector<double> > table_stableTr(90,std::vector<double>(90,-1.));
  while(std::getline(ifss2,temp_lineB))
    {
      std::stringstream stream(temp_lineB);
      stream >> nameBT >> As >> Zs >> density;
      table_stableTr[As][Zs]=density;
    }  

  
  std::ifstream ifs ( name_in.c_str() );

  std::string temp_line;
  int Af,Zf,Ap,Zp,At,Zt;
  double xs;

  std::map<int,std::map<int,std::map<int , std::map<int , std::map<int , std::map<int , double> > > > > > table_XSAll;

  std::map<std::string,std::vector<TupleElement> > BeamTarget;
  
  while(std::getline(ifs,temp_line))
    {
      
      std::stringstream stream(temp_line);
      stream >> Af >> Zf >> Ap >>Zp >> At >> Zt >> xs;

      if(table_stable[Ap][Zp]<0. || table_stableTr[At][Zt]<0)
	continue;
      //if(Zf== Zp && Af==Ap)

      double TargetDensity = table_stable[At][Zt];
      
      double count_per_cm = xs*1e-27*TargetDensity/At*6.02214129e23;
      if(not_all && count_per_cm<1e-8)
	continue;
	 
      TString nameF("");
      nameF+=Af;
      nameF+= Zf==0 ? "n" : ElName[Zf-1];
      TString nameP("");
      nameP+=Ap;
      nameP+= Zp==0 ? "n" : ElName[Zp-1];
      TString nameT(""); 
      nameT+=At;
      nameT+= Zt ==0 ? "n" : ElName[Zt-1];

      if(vecFrag.size()>0)
	{
	  int found = 0;
	  for (const std::string& FragSel : vecFrag )
	    {
	      //cout<<"test :"<<FragSel<<" "<<nameF_cp<<" cp:"<<FragSel.compare(nameF_cp)<<" "<<nameF.Contains(FragSel)<<endl;
	      if(nameF.Contains(FragSel) == true)
		{
		  found = 1;
		  break;
		}
	    }
	  if(found == 0)
	    continue;
	}
      table_XSAll[Af][Zf][Ap][Zp][At][Zt]=count_per_cm;
	
      
      TString temp_nameBeamTarget(nameP);
      temp_nameBeamTarget+="+";
      temp_nameBeamTarget+=nameT;
	 
      std::string nameBeamTarget(temp_nameBeamTarget.Data());
      std::map<std::string,std::vector<TupleElement> >::iterator it_BT = BeamTarget.find(nameBeamTarget);
      if(it_BT==BeamTarget.end())
	{
	  std::string SnameF(nameF.Data());
	  TupleElement temp_ele(SnameF,Af,Zf,1e-27*xs*TargetDensity/At*6.02214129e23);
	  std::vector<TupleElement> temp_vec(1,temp_ele);

	  BeamTarget.insert(std::pair<std::string,std::vector<TupleElement> >(nameBeamTarget,temp_vec));
	}
      else
	{
	  std::string SnameF(nameF.Data());
	  TupleElement temp_ele(SnameF,Af,Zf,1e-27*xs*TargetDensity/At*6.02214129e23);
	  it_BT->second.push_back(temp_ele);
	}
      
    }

  std::map<std::string,std::vector<std::string> > BeamTargetFragmentSister;
  std::map<std::string,std::vector<std::string> > BeamTargetFragmentSisterStrict;
  
  for(std::map<std::string,std::vector<TupleElement> >::const_iterator it_BT = BeamTarget.begin(),it_BT_end = BeamTarget.end();it_BT!=it_BT_end;++it_BT)
    {
      unsigned int size_tuple = it_BT->second.size();
      for(unsigned int id_tuple = 0;id_tuple<size_tuple;++id_tuple)
	{
	  TupleElement Tu_F(it_BT->second[id_tuple]);
	  std::string temp_key_name(it_BT->first);
	  temp_key_name+=Tu_F.name;
	  
	  for(unsigned int id_temp = 0;id_temp<size_tuple;++id_temp)
	    if(id_temp!=id_tuple)
	      {
		TupleElement Tu_temp(it_BT->second[id_temp]);
		if(Tu_F.xs<1000.*Tu_temp.xs)
		  {
		    std::map<std::string,std::vector<std::string> >::iterator it_BTF = BeamTargetFragmentSister.find(temp_key_name);
		    if(it_BTF==BeamTargetFragmentSister.end())
		      {
			std::vector<std::string> temp_vec(1,Tu_temp.name);
			BeamTargetFragmentSister.insert(std::pair<std::string,std::vector<std::string> >(temp_key_name,temp_vec));
		      }
		    else
		      {
			it_BTF->second.push_back(temp_key_name);
		      }
		  }

		if(Tu_F.xs<10.*Tu_temp.xs)
		  {
		    std::map<std::string,std::vector<std::string> >::iterator it_BTF = BeamTargetFragmentSisterStrict.find(temp_key_name);
		    if(it_BTF==BeamTargetFragmentSisterStrict.end())
		      {
			std::vector<std::string> temp_vec(1,Tu_temp.name);
			BeamTargetFragmentSisterStrict.insert(std::pair<std::string,std::vector<std::string> >(temp_key_name,temp_vec));
		      }
		    else
		      {
			it_BTF->second.push_back(temp_key_name);
		      }
		  }
	      }
	}
    }



  ofstream outfile;
  std::string outname(name_in);
  std::string ext;
  if(not_all==false)
    {
      ext = "Upated_new.dat";
    }
  else
    {
      ext = "CutSelUpated_new2.dat";
    }
  outname+=ext;
  outfile.open(outname.c_str());

  
  for(std::map<int,std::map<int,std::map<int , std::map<int , std::map<int , std::map<int , double> > > > > >::const_iterator it_Af = table_XSAll.begin(), it_Af_end = table_XSAll.end();it_Af != it_Af_end;++it_Af)
    for(std::map<int,std::map<int , std::map<int , std::map<int , std::map<int , double> > > > >::const_iterator it_Zf = it_Af->second.begin(), it_Zf_end = it_Af->second.end();it_Zf != it_Zf_end;++it_Zf)
      for(std::map<int , std::map<int , std::map<int , std::map<int , double> > > >::const_iterator it_Ap = it_Zf->second.begin(), it_Ap_end = it_Zf->second.end();it_Ap != it_Ap_end;++it_Ap)
	for(std::map<int , std::map<int , std::map<int , double> > >::const_iterator it_Zp = it_Ap->second.begin(), it_Zp_end = it_Ap->second.end();it_Zp != it_Zp_end;++it_Zp)
	  for(std::map<int , std::map<int , double> >::const_iterator it_At = it_Zp->second.begin(), it_At_end = it_Zp->second.end();it_At != it_At_end;++it_At)
	    for(std::map<int , double>::const_iterator it_Zt = it_At->second.begin(), it_Zt_end = it_At->second.end();it_Zt != it_Zt_end;++it_Zt)
	      {
		TString nameF("");
		nameF+=it_Af->first;
		nameF+= it_Zf->first ==0 ? "n" : ElName[it_Zf->first-1];
		TString nameP("");
		nameP+=it_Ap->first;
		nameP+= it_Zp->first ==0 ? "n" : ElName[it_Zp->first-1];
		TString nameT(""); 
		nameT+=it_At->first;
		nameT+= it_Zt->first ==0 ? "n" : ElName[it_Zt->first-1];

		double tempXs = it_Zt->second;
		
		double TargetDensity = table_stable[it_At->first][it_Zt->first];
		

		TString temp_nameBeamTarget(nameP);
		temp_nameBeamTarget+="+";
		temp_nameBeamTarget+=nameT;
      
		std::string nameBeamTarget(temp_nameBeamTarget.Data());
		temp_nameBeamTarget+=nameF;
		std::string nameBeamTargetFragment(temp_nameBeamTarget.Data());

		std::map<std::string,std::vector<TupleElement> >::const_iterator it_BeamTarget = BeamTarget.find(nameBeamTarget);

		std::map<std::string,std::vector<std::string> >::const_iterator it_BeamTargetFragment = BeamTargetFragmentSister.find(nameBeamTargetFragment);
		std::map<std::string,std::vector<std::string> >::const_iterator it_BeamTargetFragmentStrict = BeamTargetFragmentSisterStrict.find(nameBeamTargetFragment);
		
		int size_sister_all = it_BeamTarget!=BeamTarget.end() ? it_BeamTarget->second.size() : 0;
		int size_sister_strict1000 = it_BeamTargetFragment!=BeamTargetFragmentSister.end() ? it_BeamTargetFragment->second.size() : 0;
		int size_sister_strict10 = it_BeamTargetFragmentStrict!=BeamTargetFragmentSisterStrict.end() ? it_BeamTargetFragmentStrict->second.size() : 0;
		
		assert(it_BeamTarget != BeamTarget.end());

		outfile << it_Af->first << " "<<it_Zf->first <<  " "<<it_Ap->first << " "<< it_Zp->first << " "<< it_At->first << " "<< it_Zt->first <<" "<< TargetDensity <<" "<< tempXs*it_At->first/6.02214129e23/TargetDensity*1e27 << " "<< tempXs << " "<< size_sister_all <<" "<< size_sister_strict1000<< " "<< size_sister_strict10 << " "<<  it_Zt->first* it_Zp->first << endl ;

	      }

  outfile.close();
  return 0;
}

int UpdateFileWithFrag(std::string nameIN)
{
  //std::vector<std::string> vecF = {"18N","11C","19F","13N","21Ne","13C","14O","17F","14N","17Ne","19O","13O","18O","12C","15N","18Ne","19Ne","15C","12N"};
  //std::vector<std::string> vecF = {"12C","15C","15N","12N","13O","18O","18Ne","19Ne"};
  //std::vector<std::string> vecF = {"10C","18N","11C","13N","21Ne","13C","14O","17F","17Ne","19O"};
  std::vector<std::string> vecF = {"11B","12B","9B","10Be","11Be","12Be","14C","16C","9C","16F","18F","21F","16N","17N","22Ne","23Ne","24Ne","25Ne","15O","16O"};
  
  return UpdateFile(nameIN,"not",vecF);
}



int CorrectProductionRate(const std::string& name_database,const std::list<std::string>& list_name_in)
{


  char digit[] = "0123456789";

  std::ifstream ifs ( name_database.c_str() );

  std::string temp_line;
  int Af,Zf,Ap,Zp,At,Zt;
  double xs;
  double prodrate;
  double TargetDensity;
  std::map<std::string, std::map<std::string, std::map<std::string , double> > > table_XSAll;
  
  while(std::getline(ifs,temp_line))
    {
      
      std::stringstream stream(temp_line);

      stream >> Af >> Zf >> Ap >>Zp >> At >> Zt >> TargetDensity >> xs >> prodrate;

      TString nameF("A");
      nameF+=Af;
      nameF+="Z";
      nameF+=Zf;
      TString nameP("A");
      nameP+=Ap;
      nameP+="Z";
      nameP+=Zp;
      TString nameT("A"); 
      nameT+=At;
      nameT+="Z";
      nameT+=Zt;
      
      std::string SnameF(nameF.Data());
      std::string SnameP(nameP.Data());
      std::string SnameT(nameT.Data());

      table_XSAll[SnameF][SnameP][SnameT]=prodrate;
    }

  ifs.close();


  for(std::list<std::string>::const_iterator it_file = list_name_in.begin(),it_file_end=list_name_in.end();it_file!=it_file_end;++it_file)
    {
      
      std::string basename,dirname,namewithoutext;
      std::string::size_type p = it_file->rfind('/');
      if (p == std::string::npos) 
	{
	  std::cerr << "s contains no forward slashes" << std::endl;
	} 
      else 
	{
	  std::string s1(*it_file, 0, p);
	  std::string s2(*it_file, p + 1);
	  basename = s2;
	  dirname  = s1;
	}
 
      p = it_file->rfind('.');
      if (p == std::string::npos) 
	{
	  std::cerr << "s contains no forward dot" << std::endl;
	} 
      else 
	{
	  std::string s1(*it_file, 0, p);
	  namewithoutext = s1;
	}
  
      //cout<<" basename :"<<basename<<" "<<namewithoutext<<endl;
      stringstream stream_name(basename);
      std::string word;
      std::vector<std::string> Elements;

      while( std::getline(stream_name, word, '_') )
	{
	  if(word != "Frag" && word != "result" && word != "para.dat" && word != "result.dat")
	    {
	      //cout << word << " ";
	      //Elements.push_back(word);
	      std::string word2(word);
	  
	      for(unsigned int i = 0; i < strlen(digit); ++i)
		word.erase(std::remove(word.begin(), word.end(), digit[i]),word.end());
	  
	      string::size_type foundpos = word2.find(word);
	      if ( foundpos != string::npos )
		word2.erase(word2.begin() + foundpos, word2.begin() + foundpos + word.length());

	      //cout<< word << endl;
	      int Z = 0;
	      for(Z = 0; Z < 111 ; ++Z)
		if(word == ElName2[Z])
		  {
		    break;
		  }
	      //cout<<"A"<<word2<<"Z"<<Z<<endl;
	      TString nameTemp("A");
	      nameTemp+=word2.c_str();
	      nameTemp+="Z";
	      nameTemp+=Z;
	      Elements.push_back(std::string(nameTemp.Data()));
	    }
	}
      // for(unsigned int i=0;i<Elements.size();++i)
      // 	cout<<Elements[i]<<" ";
      // cout<<endl;


      std::ifstream ifs2 ( it_file->c_str() );

      std::getline(ifs2,temp_line);
      std::string lineComment[10];
      std::stringstream stream_init(temp_line);
      stream_init >> lineComment[0] >> lineComment[1] >> lineComment[2] >> lineComment[3] >> lineComment[4] 
		  >> lineComment[5] >> lineComment[6] >> lineComment[7] >> lineComment[8] >> lineComment[9];

      
      // cout<< table_XSAll[lineComment[3]][Elements[0]][Elements[1]] << endl;
      // cout<< table_XSAll[lineComment[5]][Elements[0]][Elements[1]] << endl;
      // cout<< table_XSAll[lineComment[7]][Elements[0]][Elements[1]] << endl;
      // cout<< table_XSAll[lineComment[9]][Elements[0]][Elements[1]] << endl;
  
      double ProdRate[4] = {    table_XSAll[lineComment[3]][Elements[0]][Elements[1]],
				table_XSAll[lineComment[5]][Elements[0]][Elements[1]],
				table_XSAll[lineComment[7]][Elements[0]][Elements[1]],  
				table_XSAll[lineComment[9]][Elements[0]][Elements[1]]  };

      
      double ProdInit = ProdRate[0];
      for(int i=0;i<4;++i)
	{
	  if(TMath::Abs(ProdRate[i])<1e-11)
	    {
	      cout<<"!>  basename :"<<basename<<" "<<namewithoutext<<" ";
	      for(unsigned int i=0;i<Elements.size();++i)
		cout<<Elements[i]<<" ";
	      cout<<endl;
	      cout<<" ["<<i<<"] "<<ProdRate[i]<<" /"<<ProdInit<<endl;
	    }
	  ProdRate[i]/=ProdInit;
	}

      // for(int i=0;i<4;++i)
      // 	cout<<ProdRate[i]<<" ";
      // cout<<endl;

      ofstream outfile;
      std::string outname(namewithoutext);
      std::string ext("_new.dat");
      outname+=ext;
      outfile.open(outname.c_str());

      for(int i=0;i<10;++i)
	{
	  //cout<<lineComment[i]<<" ";
	  outfile<<lineComment[i]<<" ";
	}     
      //cout<<endl;
      outfile<<endl;

      double Length,Thickness,Prod0,Prod1,Prod2,Prod3;
      int id0,id1,id2,id3;
      while(std::getline(ifs2,temp_line))
	{
      
	  std::stringstream stream(temp_line);
	  stream >> Length >> Thickness >> id0 >> Prod0 >> id1 >> Prod1 >> id2 >> Prod2 >> id3 >> Prod3;

	  //std::cout <<Length << " "<<Thickness << " "<<id0 << " "<<Prod0<<" "<<id1 << " "<<Prod1*ProdRate[1] << " "<<id2 << " "<<Prod2*ProdRate[2] << " "<<id3 <<" "<< Prod3*ProdRate[3] << std::endl;
	  outfile <<Length << " "<<Thickness << " "<<id0 << " "<<Prod0<<" "<<id1 << " "<<Prod1*ProdRate[1] << " "<<id2 << " "<<Prod2*ProdRate[2] << " "<<id3 <<" "<< Prod3*ProdRate[3] << std::endl;

	}
  
      outfile.close();
      ifs2.close();
      

    }
  
  return 0;
}

int CorrectOne(const std::string& name_database,const std::string& name_in)
{

  std::list<std::string> list_in;
  list_in.push_back(name_in);

  return CorrectProductionRate(name_database,list_in);

}
int CorrectAllOld(const std::string& name_database)
{

  std::list<std::string> list_in;
  list_in.push_back("../mocadi/header_files/C12_B10_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/C12_B11_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/C12_Be9_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/C12_C12_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/C12_Li7_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/C13_B10_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/C13_B11_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/C13_Be9_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/C13_C12_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/C13_Li7_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/F19_B10_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/F19_B11_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/F19_Be9_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/F19_C12_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/F19_Li7_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/N14_B10_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/N14_B11_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/N14_Be9_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/N14_C12_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/N14_Li7_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/N15_B10_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/N15_B11_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/N15_Be9_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/N15_C12_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/N15_Li7_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/Ne20_B10_Frag_C9_result_para.dat");
  list_in.push_back("../mocadi/header_files/Ne20_B11_Frag_C9_result_para.dat");
  list_in.push_back("../mocadi/header_files/Ne20_Be9_Frag_C9_result_para.dat");
  list_in.push_back("../mocadi/header_files/Ne20_C12_Frag_C9_result_para.dat");
  list_in.push_back("../mocadi/header_files/Ne20_Li7_Frag_C9_result_para.dat");
  list_in.push_back("../mocadi/header_files/O16_B10_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/O16_B11_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/O16_Be9_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/O16_C12_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/O16_Li7_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/O17_B10_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/O17_B11_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/O17_Be9_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/O17_C12_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/O17_Li7_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/O18_B10_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/O18_B11_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/O18_Be9_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/O18_C12_Frag_C9_result_para.dat"); 
  list_in.push_back("../mocadi/header_files/O18_Li7_Frag_C9_result_para.dat"); 

  return CorrectProductionRate(name_database,list_in);

}


int CorrectAllFromList(const std::string& name_database, const std::string& name_list)
{
  std::list<std::string> list_in;
  
  std::ifstream ifss(name_list.c_str());
  std::string temp_lineB,name_f;
  while(std::getline(ifss,temp_lineB))
    {
      std::stringstream stream(temp_lineB);
      stream>>name_f;
      list_in.push_back(name_f);
    }

  return CorrectProductionRate(name_database,list_in);
  
  
}


void CreateTree(const std::string& name_in, const std::string& name_out)
{
  std::ifstream ifs ( name_in.c_str() );
  
  std::string temp_line;
  Int_t Af,Zf,Ap,Zp,At,Zt;
  Float_t Rho,xs,DeXs;
  Int_t sister1,sister2,sister3,ZtZp;
  Float_t Tck;
  Float_t Survival;
  Float_t EnergyMean = 0. , EnergySigma = 0.;
  Float_t Transmit;
  Float_t Rate1, Rate2, Rate3, Rate4;
  TString nameF;
  TString nameP;
  TString nameT;
  
  TFile  *fileOut = new TFile(name_out.c_str(),"RECREATE");  
  TTree* fDataT = new TTree("EpaxDataBase","Database Epax and Mocadi");

  fDataT->Branch("Af",&Af,"Af/I");
  fDataT->Branch("Zf",&Zf,"Zf/I");
  fDataT->Branch("Ap",&Ap,"Ap/I");
  fDataT->Branch("Zp",&Zp,"Zp/I");
  fDataT->Branch("At",&At,"At/I");
  fDataT->Branch("Zt",&Zt,"Zt/I");
  fDataT->Branch("nameF",&nameF);
  fDataT->Branch("nameP",&nameP);
  fDataT->Branch("nameT",&nameT);
  fDataT->Branch("Rho",&Rho,"Rho/F");
  fDataT->Branch("xs",&xs,"xs/F");
  fDataT->Branch("DeXs",&DeXs,"DeXs/F");
  fDataT->Branch("Tck",&Tck,"Tck/F");
  fDataT->Branch("Survival",&Survival,"Survival/F");
  fDataT->Branch("EnergyMean",&EnergyMean,"EnergyMean/F");
  fDataT->Branch("EnergySigma",&EnergySigma,"EnergySigma/F");
  fDataT->Branch("Transmit",&Transmit,"Transmit/F");
  fDataT->Branch("Rate1",&Rate1,"Rate1/F");
  fDataT->Branch("Rate2",&Rate2,"Rate2/F");
  fDataT->Branch("Rate3",&Rate3,"Rate3/F");
  fDataT->Branch("Rate4",&Rate4,"Rate4/F");
  
  //std::getline(ifs,temp_line);
  while(std::getline(ifs,temp_line))
    {
      std::stringstream stream(temp_line);
      stream >> Af >> Zf >> Ap >>Zp >> At >> Zt >> Rho >> xs >> DeXs
	     >> sister1 >> sister2 >> sister3 >> ZtZp
	     >> Tck >> Survival >> EnergyMean >> EnergySigma;

      if(TMath::Abs(EnergyMean)<1e-1 || TMath::Abs(EnergySigma)<1e-1)
	continue;

      //TString nameF;
      nameF = Zf==0 ? "n" : ElName[Zf-1];
      //TString nameP;
      nameP = Zp==0 ? "n" : ElName[Zp-1];
      //TString nameT;
      nameT= Zt ==0 ? "n" : ElName[Zt-1];
    
      nameF+=Af;
      nameP+=Ap;
      nameT+=At;
      
      stream >> Transmit >> Rate1 >> Rate2 >> Rate3 >> Rate4;

      fDataT->Fill();

    }
  
  fileOut->Write();
  
  
}


struct tuple_graph 
{
public :
  std::string BT;

  double TargetGramCm;
  double SurvivalRate;
  double Transmittance;
  double EnergyMean;
  double EnergySigma;
  double ProdFrag;
  double ProdPara1;
  double ProdPara2;
  double ProdPara3;

    
//  tuple_graph()
//  {}
// tuple_graph(const std::string& Name, double TgCM, double prodF, double prodP1, double prodP2, double prodP3):
//   BT(Name),TargetCm(1,TCm),TargetGramCm(1,TgCM),ProdFrag(1,prodF),ProdPara1(1,prodP1),ProdPara2(1,prodP2),ProdPara3(1,prodP3),ProdStage(1,stage)
// {
  
// }
//~tuple_graph() {}
// tuple_graph(const tuple_graph& t):BT(t.BT),TargetCm(t.TargetCm),TargetGramCm(t.TargetGramCm),ProdFrag(t.ProdFrag),ProdPara1(t.ProdPara1),ProdPara2(t.ProdPara1),ProdPara3(t.ProdPara3),ProdStage(t.ProdStage)
// {}
// void AddValue(double TCm, double TgCM, double prodF, double prodP1, double prodP2, double prodP3,const std::vector<double>& stage)
// {
//   TargetCm.push_back(TCm);
//   TargetGramCm.push_back(TgCM);
//   ProdFrag.push_back(prodF);
//   ProdPara1.push_back(prodP1);
//   ProdPara2.push_back(prodP2);
//   ProdPara3.push_back(prodP3);
//   ProdStage.push_back(stage);
// }
};


void Plots(int Z, int A, const std::string& name_in, double Eth = 2000.)
{

  TGraphErrors* GraphEnergy = new TGraphErrors();

  std::ifstream ifs ( name_in.c_str() );
  
  std::string temp_line;
  Int_t Af,Zf,Ap,Zp,At,Zt;
  Float_t Rho,xs,DeXs;
  Int_t sister1,sister2,sister3,ZtZp;
  Float_t Tck;
  Float_t Survival;
  Float_t EnergyMean = 0. , EnergySigma = 0.;
  Float_t Transmit;
  Float_t Rate1, Rate2, Rate3, Rate4;
  TString nameF;
  TString nameP;
  TString nameT;

  std::set<std::string> allFrag;

  std::map<std::string, std::map<int, tuple_graph> > BeamTarget;
  
  while(std::getline(ifs,temp_line))
    {
      std::stringstream stream(temp_line);
      stream >> Af >> Zf;

      //TString nameF;
      nameF = Zf==0 ? "n" : ElName[Zf-1];
      nameF+=Af;

      std::string SnameF(nameF.Data());
      allFrag.insert(SnameF);
      
      
      stream >> Ap >>Zp >> At >> Zt >> Rho >> xs >> DeXs
	     >> sister1 >> sister2 >> sister3 >> ZtZp
	     >> Tck >> Survival >> EnergyMean >> EnergySigma;

      GraphEnergy->SetPoint(GraphEnergy->GetN(),GraphEnergy->GetN(),EnergyMean);
      GraphEnergy->SetPointError(GraphEnergy->GetN()-1,0,EnergySigma);
      if(Af != A || Zf != Z)
	continue;
      
      if(Zt == 1) 
	continue;

      if(Zt != 4 && Zt != 5 && Zt != 6)
	continue;
      
      if(TMath::Abs(EnergyMean)<1e-1 || TMath::Abs(EnergySigma)<1e-1)
	continue;

      stream >> Transmit >> Rate1 >> Rate2 >> Rate3 >> Rate4;

      //TString nameP;
      nameP = Zp==0 ? "n" : ElName[Zp-1];
      //TString nameT;
      nameT= Zt ==0 ? "n" : ElName[Zt-1];
    
      nameP+=Ap;
      nameT+=At;
      
      std::string SnameP(nameP.Data());
      std::string SnameT(nameT.Data());
      std::string SnameAll(SnameF);
      SnameAll+=SnameP;
      SnameAll+=SnameT;
      std::string SnameBT(SnameP);
      SnameBT+="+";
      SnameBT+=SnameT;

      
      
      int Thick = static_cast<int>(Tck);
      
      auto it_FBT = BeamTarget.find(SnameAll);
      if(it_FBT==BeamTarget.end())
	{
	  std::map<int, tuple_graph> tempMap;
	  
	  tuple_graph temp = {SnameBT, Rho, Survival, Transmit, EnergyMean, EnergySigma, Rate1, Rate2, Rate3, Rate4};
	  tempMap.insert(std::pair<int,tuple_graph>(Thick,temp));
	  BeamTarget.insert(std::pair<std::string,std:: map<int, tuple_graph> >(SnameAll,tempMap));
	}
      else
	{
	  auto it_BT = it_FBT->second.find(Thick);
	  if(it_BT == it_FBT->second.end())
	    {
	      tuple_graph temp = {SnameBT, Rho, Survival, Transmit, EnergyMean, EnergySigma, Rate1, Rate2, Rate3, Rate4};
	      it_FBT->second.insert(std::pair<int,tuple_graph>(Thick,temp));
	    }
	  else
	    {
	      cout<<"!> Already inserted line ! "<<SnameAll<<" "<<Tck<<endl;
	    }
	}
    }


  for(auto it_F = allFrag.begin(), it_F_end = allFrag.end();it_F!=it_F_end;++it_F)
    {
      cout<<*it_F<<endl;
    }

  // for (auto it_FBT = BeamTarget.begin(), it_FBT_end = BeamTarget.end(); it_FBT != it_FBT_end; ++it_FBT)
  //   {
  //     cout<<it_FBT->first<<"\n";
  //     for(auto it_BT = it_FBT->second.begin(), it_BT_end = it_FBT->second.end(); it_BT != it_BT_end; ++it_BT)
  // 	{
  // 	  cout<<"--: Tck "<<it_BT->first<<" "<<it_BT->second.BT<<" "<<it_BT->second.TargetGramCm<<" "<<it_BT->second.SurvivalRate<<" "<<it_BT->second.EnergyMean<<" "<<it_BT->second.ProdFrag<<"\n";
	 
  // 	}
  //     cout<<" ----- "<<endl;
  //   }
  
  TGraph* GraphAll[BeamTarget.size()] ;
  TGraph* GraphAllRatio[BeamTarget.size()];

  //TGraph* GraphAllRatio[BeamTarget.size()];

  TMultiGraph* GraphM = new TMultiGraph();
  TMultiGraph* GraphRatioM = new TMultiGraph();

  int i = 0;
  for(auto it_FBT = BeamTarget.begin(), it_FBT_end = BeamTarget.end();it_FBT!= it_FBT_end;++it_FBT)
    {
      //cout<<"Status :"<<it_FBT->first<<" | "<<it_FBT->second.BT<<" "<<it_FBT->second.TargetCm.size()<<endl;

      GraphAll[i] = new TGraph();
      GraphAllRatio[i] = new TGraph();

      auto it_BT = it_FBT->second.begin();
      
      GraphAll[i]->SetNameTitle(it_BT->second.BT.c_str(),it_BT->second.BT.c_str());
      GraphAllRatio[i]->SetNameTitle(it_BT->second.BT.c_str(),it_BT->second.BT.c_str());

      int k = 0;
      for(auto it_BT_end = it_FBT->second.end();it_BT!=it_BT_end;++it_BT)
	{
	  
	  if(it_BT->second.EnergyMean>Eth  && it_BT->first>0)
	    {
	      GraphAll[i]->SetPoint(k,it_BT->second.TargetGramCm*it_BT->first,1e10*it_BT->second.ProdFrag*(1.-it_BT->second.SurvivalRate));

	      double temp_ratio = TMath::Abs(it_BT->second.ProdPara1+it_BT->second.ProdPara3)>1e-6 ? 1e10*(it_BT->second.ProdFrag-(it_BT->second.ProdPara1+it_BT->second.ProdPara3))*(1.-it_BT->second.SurvivalRate) : 1e10*it_BT->second.ProdFrag*(1.-it_BT->second.SurvivalRate);

	      GraphAllRatio[i]->SetPoint(k,it_BT->second.TargetGramCm*it_BT->first,temp_ratio);
	      ++k;
	    }
	}
      auto it_BT2 = it_FBT->second.rbegin();
      if(1e10*it_BT2->second.ProdFrag*(1.-it_BT2->second.SurvivalRate)>0.5e6)
	{	  
	  GraphAll[i]->SetLineColor(i%9+1);
	  GraphAllRatio[i]->SetLineColor(i%9+1);
	}
      ++i;
      
      
      //GraphM->Add(GraphAll[i]);
      //GraphRatioM->Add(GraphAllRatio[i]);
    }

  int NN = 0;
  double MaxBeam = -1;
  double MinBeam = 10000;
  double MeanMaxBeam = -1;
  for(unsigned int k = 0;k<BeamTarget.size();++k)
    {
      if(0==GraphAll[k])
	cout<<"E> GrapAll ["<<k<<"] NULL !"<<endl;
      else
	{
	  //cout<<"id#"<<k<<endl; GraphAll[k]->Print();
	  int tempN = GraphAll[k]->GetN();
	  double tempX,tempY;
	  GraphAll[k]->GetPoint(tempN-1,tempX,tempY);
	  double tempY2;
	  GraphAll[k]->GetPoint(tempN/2,tempX,tempY2);
	  
	  if(MaxBeam<tempY)
	    MaxBeam = tempY;
	  if(MinBeam>tempY)
	    MinBeam = tempY;
	  if(MeanMaxBeam<tempY2)
	    MeanMaxBeam = tempY;
	}
    }

  int Nx = 0;
  for(unsigned int k = 0;k<BeamTarget.size();++k)
    {
      if(0==GraphAll[k])
	cout<<"E> GrapAll ["<<k<<"] NULL !"<<endl;
      else
	{
	  //cout<<"id#"<<k<<endl; GraphAll[k]->Print();
	  int tempN = GraphAll[k]->GetN();
	  double tempX,tempY;
	  GraphAll[k]->GetPoint(tempN-1,tempX,tempY);
	  if(tempY>MaxBeam*0.80)
	    {
	      GraphAll[k]->SetMarkerStyle(1);
	      GraphM->Add(GraphAll[k],"C");
	      ++NN;
	    }
	  else if(tempY<MinBeam*1.1)
	    {
	      double rand = gRandom->Uniform(0,1);
	      int Ntest = 10*rand;
	      if(Nx%10==Ntest)
		{
		  GraphAll[k]->SetMarkerStyle(1);
		  GraphM->Add(GraphAll[k],"C");
		  ++NN;
		}
	    }
	  else if(tempY>MeanMaxBeam*0.8)
	    {
	      double rand = gRandom->Uniform(0,1);
	      int Ntest = 10*rand;
	      // if(Nx%10==Ntest)
		{
		  GraphAll[k]->SetMarkerStyle(1);
		  GraphM->Add(GraphAll[k],"C");
		  ++NN;
		}

	    }
	  else
	    {
	      ++Nx;
	      double rand = gRandom->Uniform(0,1);
	      int Ntest = 10*rand;
	      if(Nx%10==Ntest)
		{
		  GraphAll[k]->SetMarkerStyle(1);
		  GraphM->Add(GraphAll[k],"C");
		  ++NN;
		}
	    }
	    
	}
      if(0==GraphAllRatio[k])
	cout<<"E> GrapAllRatio ["<<k<<"] NULL !"<<endl;
      else
	{
	  GraphRatioM->Add(GraphAllRatio[k]);
	}
    }
  cout<<" N:"<<BeamTarget.size()<<" "<<NN<<endl;
  
  TCanvas* c1 = new TCanvas("c1","c1",500,500);
  c1->cd();
  //GraphAll[0]->Draw("A*L");
  GraphM->Draw("a*c");
  TCanvas* c2 = new TCanvas("c2","c2",500,500);
  c2->cd();
  //GraphAllRatio[0]->Draw("A*L");
  GraphRatioM->Draw("a*c");

  TCanvas* c3 = new TCanvas("c3","c3",500,500);
  c3->cd();
  GraphEnergy->Draw("ap");
  
}
