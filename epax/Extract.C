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

void ExtractAll(std::string name_in,bool updated=false)
{

  std::ifstream ifss("");


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
	  h_Fxs->Fill(N,Z,table_XS[Z][A]);
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
  h_Fxs->Draw("colz text");
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

class tuple_element 
{
public :
  std::string name;
  int A;
  int Z;
  double xs;
  
  tuple_element():name(""),A(0),Z(0),xs(-1)
  {}
  tuple_element(std::string N,int Af, int Zf, double XS):name(N),A(Af),Z(Zf),xs(XS)
  {}
  ~tuple_element() {}
  tuple_element(const tuple_element& t):name(t.name),A(t.A),Z(t.Z),xs(t.xs)
  {}
};

int UpdateFile(std::string name_in,std::string option="all")
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

  std::ifstream ifs ( name_in.c_str() );

  std::string temp_line;
  int Af,Zf,Ap,Zp,At,Zt;
  double xs;

  std::map<int,std::map<int,std::map<int , std::map<int , std::map<int , std::map<int , double> > > > > > table_XSAll;

  std::map<std::string,std::vector<tuple_element> > BeamTarget;
  
  while(std::getline(ifs,temp_line))
    {
      
      std::stringstream stream(temp_line);
      stream >> Af >> Zf >> Ap >>Zp >> At >> Zt >> xs;

      if(table_stable[Ap][Zp]<0. || table_stable[At][Zt]<0)
	continue;
      //if(Zf== Zp && Af==Ap)

      double TargetDensity = table_stable[At][Zt];
      
      double count_per_cm = xs*1e-27*TargetDensity/At*6.02214129e23;
      if(not_all && count_per_cm<1e-8)
	continue;
	 
      table_XSAll[Af][Zf][Ap][Zp][At][Zt]=count_per_cm;
	
      TString nameF("");
      nameF+=Af;
      nameF+= Zf==0 ? "n" : ElName[Zf-1];
      TString nameP("");
      nameP+=Ap;
      nameP+= Zp==0 ? "n" : ElName[Zp-1];
      TString nameT(""); 
      nameT+=At;
      nameT+= Zt ==0 ? "n" : ElName[Zt-1];
      
      TString temp_nameBeamTarget(nameP);
      temp_nameBeamTarget+="+";
      temp_nameBeamTarget+=nameT;
	 
      std::string nameBeamTarget(temp_nameBeamTarget.Data());
      std::map<std::string,std::vector<tuple_element> >::iterator it_BT = BeamTarget.find(nameBeamTarget);
      if(it_BT==BeamTarget.end())
	{
	  std::string SnameF(nameF.Data());
	  tuple_element temp_ele(SnameF,Af,Zf,1e-27*xs*TargetDensity/At*6.02214129e23);
	  std::vector<tuple_element> temp_vec(1,temp_ele);

	  BeamTarget.insert(std::pair<std::string,std::vector<tuple_element> >(nameBeamTarget,temp_vec));
	}
      else
	{
	  std::string SnameF(nameF.Data());
	  tuple_element temp_ele(SnameF,Af,Zf,1e-27*xs*TargetDensity/At*6.02214129e23);
	  it_BT->second.push_back(temp_ele);
	}
      
    }

  std::map<std::string,std::vector<std::string> > BeamTargetFragmentSister;
  std::map<std::string,std::vector<std::string> > BeamTargetFragmentSisterStrict;
  
  for(std::map<std::string,std::vector<tuple_element> >::const_iterator it_BT = BeamTarget.begin(),it_BT_end = BeamTarget.end();it_BT!=it_BT_end;++it_BT)
    {
      unsigned int size_tuple = it_BT->second.size();
      for(unsigned int id_tuple = 0;id_tuple<size_tuple;++id_tuple)
	{
	  tuple_element Tu_F(it_BT->second[id_tuple]);
	  std::string temp_key_name(it_BT->first);
	  temp_key_name+=Tu_F.name;
	  
	  for(unsigned int id_temp = 0;id_temp<size_tuple;++id_temp)
	    if(id_temp!=id_tuple)
	      {
		tuple_element Tu_temp(it_BT->second[id_temp]);
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
      ext = "Upated.dat";
    }
  else
    {
      ext = "CutSelUpated.dat";
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

		std::map<std::string,std::vector<tuple_element> >::const_iterator it_BeamTarget = BeamTarget.find(nameBeamTarget);

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
  
      cout<<" basename :"<<basename<<" "<<namewithoutext<<endl;
      stringstream stream_name(basename);
      std::string word;
      std::vector<std::string> Elements;

      while( std::getline(stream_name, word, '_') )
	{
	  if(word != "Frag" && word != "result" && word != "para.dat" && word != "result.dat")
	    {
	      cout << word << " ";
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
	      cout<<"A"<<word2<<"Z"<<Z<<endl;
	      TString nameTemp("A");
	      nameTemp+=word2.c_str();
	      nameTemp+="Z";
	      nameTemp+=Z;
	      Elements.push_back(std::string(nameTemp.Data()));
	    }
	}
      for(unsigned int i=0;i<Elements.size();++i)
	cout<<Elements[i]<<" ";
      cout<<endl;


      std::ifstream ifs2 ( it_file->c_str() );

      std::getline(ifs2,temp_line);
      std::string lineComment[10];
      std::stringstream stream_init(temp_line);
      stream_init >> lineComment[0] >> lineComment[1] >> lineComment[2] >> lineComment[3] >> lineComment[4] 
		  >> lineComment[5] >> lineComment[6] >> lineComment[7] >> lineComment[8] >> lineComment[9];

 
      cout<< table_XSAll[lineComment[3]][Elements[0]][Elements[1]] << endl;
      cout<< table_XSAll[lineComment[5]][Elements[0]][Elements[1]] << endl;
      cout<< table_XSAll[lineComment[7]][Elements[0]][Elements[1]] << endl;
      cout<< table_XSAll[lineComment[9]][Elements[0]][Elements[1]] << endl;
  
      double ProdRate[4] = {    table_XSAll[lineComment[3]][Elements[0]][Elements[1]],
				table_XSAll[lineComment[5]][Elements[0]][Elements[1]],
				table_XSAll[lineComment[7]][Elements[0]][Elements[1]],  
				table_XSAll[lineComment[9]][Elements[0]][Elements[1]]  };
  
      double ProdInit = ProdRate[0];
      for(int i=0;i<4;++i)
	ProdRate[i]/=ProdInit;

      for(int i=0;i<4;++i)
	cout<<ProdRate[i]<<" ";
      cout<<endl;

      ofstream outfile;
      std::string outname(namewithoutext);
      std::string ext("_new.dat");
      outname+=ext;
      outfile.open(outname.c_str());

      for(int i=0;i<10;++i)
	{
	  cout<<lineComment[i]<<" ";
	  outfile<<lineComment[i]<<" ";
	}     
      cout<<endl;
      outfile<<endl;

      double Length,Thickness,Prod0,Prod1,Prod2,Prod3;
      int id0,id1,id2,id3;
      while(std::getline(ifs2,temp_line))
	{
      
	  std::stringstream stream(temp_line);
	  stream >> Length >> Thickness >> id0 >> Prod0 >> id1 >> Prod1 >> id2 >> Prod2 >> id3 >> Prod3;

	  std::cout <<Length << " "<<Thickness << " "<<id0 << " "<<Prod0<<" "<<id1 << " "<<Prod1*ProdRate[1] << " "<<id2 << " "<<Prod2*ProdRate[2] << " "<<id3 <<" "<< Prod3*ProdRate[3] << std::endl;
	  outfile <<Length << " "<<Thickness << " "<<id0 << " "<<Prod0<<" "<<id1 << " "<<Prod1*ProdRate[1] << " "<<id2 << " "<<Prod2*ProdRate[2] << " "<<id3 <<" "<< Prod3*ProdRate[3] << std::endl;

	}
  
      outfile.close();
      ifs2.close();
      

    }
  
  return 0;
}

int CorrectAll(const std::string& name_database)
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
