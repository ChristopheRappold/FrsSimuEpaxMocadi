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





int SimplifyData(const std::string& name_in,const std::string& name_dir)
{

  ofstream outfile;
  std::string outname(name_in);
  std::string ext("_new.dat");
  outname+=ext;
  outfile.open(outname.c_str());
 
  std::ifstream ifs ( name_in.c_str() );

  std::string temp_line;
  int Af,Zf,Ap,Zp,At,Zt;
  double xs;
  double prodrate;
  double TargetDensity;
  int Sister1,Sister2,Sister3,Z2;
  std::string name_file1;
  std::string name_file2;
  //std::string name_dir("../mocadi/");

  while(std::getline(ifs,temp_line))
    {
      
      std::stringstream stream(temp_line);
      stream >> Af >> Zf >> Ap >>Zp >> At >> Zt >> TargetDensity >> xs >> prodrate >> Sister1 >> Sister2 >> Sister3 >> Z2 >> name_file1 >> name_file2; 

      cout<<"process file "<<name_file1<<" "<<name_file2<<endl;
      std::string temp_name_file1(name_dir);
      temp_name_file1+= name_file1;
      std::string temp_name_file2(name_dir);
      temp_name_file2+= name_file2;
      
      std::ifstream temp_if1 ( temp_name_file1.c_str() );
      
      std::vector< double > Dtarget;
      std::vector< double > DtargetDensity;
      std::vector< std::vector<double > > Dres;

 
      while(std::getline(temp_if1,temp_line))
	{
	  double targetCm, targetCmDensity;
	  int Nexp;
	  std::stringstream stream2(temp_line);
	  stream2 >> targetCm >> targetCmDensity >> Nexp;
	  Dtarget.push_back(targetCm);
	  DtargetDensity.push_back(targetCmDensity);
	  std::vector<double> temp_res;
	  for(int k = 0; k<Nexp;++k)
	    {
	      double temp_pourcent, temp_prod;
	      stream2 >> temp_pourcent >> temp_prod;
	      temp_res.push_back(temp_pourcent);
	      temp_res.push_back(temp_prod);
	    }
	  Dres.push_back(temp_res);
	}

      temp_if1.close();

      std::ifstream temp_if2 ( temp_name_file2.c_str() );
      
      std::vector< double > Dtarget2;
      std::vector< double > DtargetDensity2;
      std::vector< std::vector<double > > Dres2;
      std::getline(temp_if2,temp_line);
	
      while(std::getline(temp_if2,temp_line))
	{
	  std::vector<double> temp_res;

	  double targetCm, targetCmDensity;
	  int N1,N2,N3,N4;
	  double res1,res2,res3,res4;
	  std::stringstream stream3(temp_line);
	  stream3 >> targetCm >> targetCmDensity >> N1 >> res1 >> N2 >> res2 >> N3 >> res3 >> N4 >> res4;
	  Dtarget2.push_back(targetCm);
	  DtargetDensity2.push_back(targetCmDensity);

	  temp_res.push_back(N1);
	  temp_res.push_back(res1);
	  temp_res.push_back(N2);
	  temp_res.push_back(res2);
	  temp_res.push_back(N3);
	  temp_res.push_back(res3);
	  temp_res.push_back(N4);
	  temp_res.push_back(res4);
	  Dres2.push_back(temp_res);
	}

      temp_if2.close();

      for(unsigned int i=0;i<Dtarget.size();++i)
	{
	  outfile << Af<<" " << Zf<<" " << Ap<<" " <<Zp<<" " << At<<" " << Zt<<" " << TargetDensity<<" " << xs<<" " << prodrate<<" " << Sister1<<" " << Sister2<<" " << Sister3<<" " << Z2<<" ";
	  outfile << Dtarget[i]<<" " << DtargetDensity[i]<<" "<<Dres[i].size() ;
	  for(unsigned int k=0;k<Dres[i].size();++k)
	    {
	      outfile<<" "<<Dres[i][k];
	    }
	  outfile<<" "<<Dres2[i].size();
	  for(unsigned int k=0;k<Dres2[i].size();++k)
	    {
	      outfile<<" "<<Dres2[i][k];
	    }
	  outfile<<endl;
	}
    }

  outfile.close();


  return 0;



}
