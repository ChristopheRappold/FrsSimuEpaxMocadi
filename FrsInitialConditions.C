#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>

#include "TMath.h"
#include "TString.h"
#include "TH1F.h"

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




int Solve(const std::vector<double>& p,std::vector<double>& x_out)
{
  double delta = p[1]*p[1]- 4.*p[0]*p[2];
  if(delta>0)
    {
      double temp = -p[1]/(2.*p[0]);
      double temp2 = TMath::Sqrt(delta)/(2.*p[0]);
      
      double x_t1[2] = {temp+temp2,temp-temp2};
      x_out.resize(2);
     if(x_t1[0]<x_t1[1])
	{x_out[0]=x_t1[0];x_out[1]=x_t1[1];}
      else
	{x_out[0]=x_t1[1];x_out[1]=x_t1[0];}
      return 1;
    }
  else
    return 0;
}

void SIS_Beam(std::string name_in,double EkI=2.,double BrhoI = 18.)
{




  std::ifstream ifss("epax/ListBeam.dat");


  std::ifstream ifs ( name_in.c_str() );
  std::vector<std::vector<bool> > table_stable(50,std::vector<bool>(50,false));
  std::string temp_lineB,nameBT;
  int As,Zs;
  while(std::getline(ifss,temp_lineB))
    {
      std::stringstream stream(temp_lineB);
      stream >> nameBT >> As >> Zs;
      table_stable[As][Zs]=true;
    }  

  std::string temp_line;
  int A,Z,Q;
  std::string El;
  int other0,other1;
  double MoverQ;

  std::getline(ifs,temp_line);
  bool new_type(false);
  double u = 0.931494061;
  std::map<int,std::map<int,std::map<int,double> > > table_beam;
  int beam_passed=0;
  int beam_all=0;


  
  TH1F* h_Ek = new TH1F("Energy_K","Energy_K",200,0,200);
  TH1F* h_Brho = new TH1F("Brho","Brho",200,0,200);
  

  while(std::getline(ifs,temp_line))
    {
      ++beam_all;
      std::stringstream stream(temp_line);
      if(!new_type)
	stream >> A >> Z >> Q >> El >> other0 >> MoverQ;
      else
	stream >> A >> Z >> Q >> other1 >> El >> other0 >> MoverQ;
      
      table_beam[A][Z][Q] = MoverQ;

      if(A==257 && Z==99 && Q==99)
	new_type=true;

      
      
      if(Z>15 || table_stable[A][Z]==false)
	continue;
      double MassU = MoverQ*Q;
      double Mass = MassU*u;
      double Energy = EkI*A + MassU*u;
      double Du = Energy / MassU ;





      double p = TMath::Sqrt((EkI*A+Mass)*(EkI*A+Mass) - Mass*Mass);
      double Brho = 3.10715497*p/Q;

      double params[3] = {-A*A,-A*Mass,(BrhoI*Q/3.10715497)*(BrhoI*Q/3.10715497)};
      std::vector<double> Paras,Ek_out;
      Paras.assign(params,params+3);
      
      if(Solve(Paras,Ek_out))
	{
	  //std::cout<<"A"<<A<<"Z"<<Z<<" (with charge "<<Q<<")"<<MoverQ<<" "<<Brho<<std::endl;
	  double Ek2 = Ek_out[0]>0 ? Ek_out[0] : Ek_out[1];
	  TString name_temp(ElName[Z-1].c_str());
	  name_temp+=A;
	  name_temp+="_Q";
	  name_temp+=Q;
	  
	  if(Ek2>2.)
	    {
	      ++beam_passed;
	      std::cout<<name_temp<<" A"<<A<<"Z"<<Z<<" (with charge "<<Q<<")"<<MoverQ<<" "<< Ek2<<std::endl;
	      //  }
	  //if(Z<5)
	  // {
	      h_Ek->Fill(name_temp,Ek2);
	    }
	}

      if(Bhro<=18.)
	{
	  h_Brho->Fill(name_temp,Bhro);
	}
    }

  std::cout<<beam_passed<<" / "<<beam_all<<endl;
  h_Ek->Draw();
}
