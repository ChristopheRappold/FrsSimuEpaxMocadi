#ifndef FRSANA_H
#define FRSANA_H

#include <vector>
#include <list>
#include <string>
#include <iostream>
#include <algorithm>

#include "TMath.h"
#include "TRandom3.h"
#include "TGeoManager.h"
#include "TGeoElement.h"

class Compute
{
public:
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

	//std::cout<<" Compute MeanSigma:"<<std::endl;
	// for(auto item : vec)
	// 	std::cout<<" "<<item;
	// std::cout<<std::endl;
	//std::cout<<" Mean "<<res[0]<<" Sigma "<<res[1]<<std::endl;

	return;
      }
  }

  void computeQuantile(const std::vector<double>& vec, double* res)
  {
    Double_t ProbQuantile[3] = {0.25,0.50,0.75};
    Double_t ProbRes[3] = {0,0,0};
    double Vec[ vec.size()] ;
    std::copy(vec.begin(), vec.end(), Vec);
    TMath::Quantiles(vec.size(), 3, Vec, ProbRes, ProbQuantile, false);
  
    res[0] = ProbRes[0];
    res[1] = ProbRes[1];
    res[2] = ProbRes[2];
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
    double in2 = ComputeWeigth(in[0],in[1],NtypeComp);
    if(TMath::Abs(in2)<1e-4)
      return -2;
    if(type=="A B AB")
      {
	out[0] = in[0];
	out[1] = in[1];
	out[2] = in2;//1-in[0]-in[1];
      }
    else if(type=="B A AB")
      {
	out[0] = in[1];
	out[1] = in[0];
	out[2] = in2;//1-in[0]-in[1];
      }
    else if(type=="A AB B")
      {
	out[0] = in[0];
	out[1] = in2;//1-in[0]-in[1];
	out[2] = in[1];
      }
    else if(type=="B AB A")
      {
	out[0] = in[1];
	out[1] = in2;//1-in[0]-in[1];
	out[2] = in[0];
      }
    else if(type=="AB A B")
      {
	out[0] = in2;//1-in[0]-in[1];
	out[1] = in[0];
	out[2] = in[1];
      }
    else if(type=="AB B A")
      {
	out[0] = in2;//1-in[0]-in[1];
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



  void Invariant(std::vector<double>& In, const std::vector<double*>& Range, int type, TH2F* h_hist = nullptr)
  {
    for(size_t i = 0; i<In.size();++i)
      {
	if(type == 0)
	  {
	    In[i] -= Range[i][0];
	    In[i] /= Range[i][1]-Range[i][0];
	  }
	else if(type == 1)
	  {
	    In[i] -= Range[i][0];
	    In[i] /= Range[i][1];
	  }
	else if(type == 2)
	  {
	    In[i] -= Range[i][0];
	    In[i] /= 0.5*(Range[i][2]-Range[i][1]);

	  }
	if(h_hist!=nullptr)
	  h_hist->Fill(6+i,In[i]);
      }
  }

  TGeoManager *geom; 
  TGeoElementTable *tableRN; 
  std::map<std::string,double> cache_NucleiMass;

  Compute():geom(new TGeoManager("geom","radionuclides")),tableRN(geom->GetElementTable())
  { }
  
  double Brho(int Z, int A, double Ek)
  {
    if(geom==nullptr)
       geom = new TGeoManager("geom","radionuclides");
    if(tableRN==nullptr)
      tableRN = geom->GetElementTable();
    const double u = 0.931494061;
    const double m_eminus = 5.10998909999999971e-04; // GeV
    std::string nameElem = std::to_string(Z);
    nameElem+="_";
    nameElem+=std::to_string(A);
    
    double Mass  = 0.;
    auto it_Elem = cache_NucleiMass.find(nameElem);
    if(it_Elem == cache_NucleiMass.end())
      {
	auto* TempElement = tableRN->GetElementRN(A,Z);
	double Dmass = 0.;
	if(TempElement==nullptr)
	  {
	    cout<<"E> no element ! "<<A<<" "<<Z<<" "<<TempElement<<endl;
	    if(A==23 && Z==14)
	      Dmass = 23.073*1e-3; // MeV -> GeV
	    else
	      return 999999;
	  }
	else
	  Dmass = TempElement->MassEx()*1e-3; // MeV -> GeV
	
	Mass = Dmass+A*u - Z*m_eminus;
	cache_NucleiMass.emplace(nameElem,Mass);
      }
    else
      Mass = it_Elem->second;
      
    double EkI = Ek > 50 ? Ek*1e-3 : Ek; //  MeV -> GeV  
    double p = TMath::Sqrt((EkI*A+Mass)*(EkI*A+Mass) - Mass*Mass);
    const double inv_c = 3.33564095198152044;
    double Brho = inv_c*p/Z;
    return Brho;
    
  }
  
};

class LambdaCrossSectionFunction
{
  double norm;
  double C;
  double alpha;
  double Function(double E) const
  {
    double E1 = E < 10 ? E*1e3 : E ; 
    double s = 938.272*938.272 + 938.272*938.272 + 2*938.272*(E1+938.272);
    double e = TMath::Sqrt(s) -  2.547629e3; // sqrt(s) - mp - mLambda - mK+    
    return C*e*e/(1+TMath::Sqrt(1+e/alpha))/(1+TMath::Sqrt(1+e/alpha));
  }
public:
  LambdaCrossSectionFunction(double Enorm_=-1):norm(1),C(0.02574),alpha(5.203)
  {
    if(Enorm_ > 0)
      norm = Function(Enorm_);
  }
  double operator() (double E) const
  {
    double temp = Function(E)/norm;
    return TMath::IsNaN(temp) ? 0. : temp;
  }
  double operator() (Double_t* x, Double_t* par)
  {
    return this->operator()(x[0]);
  }
};

class tuple_graph 
{
public :
  std::string BT;
  std::string F;
  std::vector<int> BT_ids;
  std::vector<double> TargetCm;
  std::vector<double> TargetGramCm;
  
  std::vector<double> ProdFrag;
  std::vector<double> ProdPara1;
  std::vector<double> ProdPara2;
  std::vector<double> ProdPara3;

  std::vector< std::vector<double> > ProdStage;
  
  std::vector<double> EnergyMean;
  std::vector<double> EnergySigma;
  std::vector<double> Survival; 

  tuple_graph()
  {}
  tuple_graph(const std::string& Name, const std::string& NameF, const std::vector<int>& ids, double TCm, double TgCM, double prodF, double prodP1, double prodP2, double prodP3,const std::vector<double>& stage, double EnergyM = 0, double EnergyS = 0, double SurvivalRate = 0):
    BT(Name),F(NameF),BT_ids(ids),TargetCm(1,TCm),TargetGramCm(1,TgCM),ProdFrag(1,prodF),ProdPara1(1,prodP1),ProdPara2(1,prodP2),ProdPara3(1,prodP3),ProdStage(1,stage),EnergyMean(1,EnergyM),EnergySigma(1,EnergyS),Survival(1,SurvivalRate)
  {}
  tuple_graph(const std::string& Name,const std::string& NameF , const std::vector<int>& ids, double TCm, double TgCM, double prodF, double prodP1, double prodP2, double prodP3, double EnergyM = 0, double EnergyS = 0, double SurvivalRate = 0):
    BT(Name),F(NameF),BT_ids(ids),TargetCm(1,TCm),TargetGramCm(1,TgCM),ProdFrag(1,prodF),ProdPara1(1,prodP1),ProdPara2(1,prodP2),ProdPara3(1,prodP3),EnergyMean(1,EnergyM),EnergySigma(1,EnergyS),Survival(1,SurvivalRate)
  {}
  ~tuple_graph() {}
  tuple_graph(const tuple_graph& t):BT(t.BT),F(t.F),BT_ids(t.BT_ids),TargetCm(t.TargetCm),TargetGramCm(t.TargetGramCm),ProdFrag(t.ProdFrag),ProdPara1(t.ProdPara1),ProdPara2(t.ProdPara1),ProdPara3(t.ProdPara3),ProdStage(t.ProdStage),EnergyMean(t.EnergyMean),EnergySigma(t.EnergySigma),Survival(t.Survival)
  {}
  void AddValue(double TCm, double TgCM, double prodF, double prodP1, double prodP2, double prodP3,const std::vector<double>& stage, double EnergyM = 0, double EnergyS = 0, double SurvivalRate = 0 )
  {
    TargetCm.push_back(TCm);
    TargetGramCm.push_back(TgCM);
    ProdFrag.push_back(prodF);
    ProdPara1.push_back(prodP1);
    ProdPara2.push_back(prodP2);
    ProdPara3.push_back(prodP3);
    ProdStage.push_back(stage);
    EnergyMean.push_back(EnergyM);
    EnergySigma.push_back(EnergyS);
    Survival.push_back(SurvivalRate);
  }
  void AddValue(double TCm, double TgCM, double prodF, double prodP1, double prodP2, double prodP3, double EnergyM = 0, double EnergyS = 0, double SurvivalRate = 0)
  {
    TargetCm.push_back(TCm);
    TargetGramCm.push_back(TgCM);
    ProdFrag.push_back(prodF);
    ProdPara1.push_back(prodP1);
    ProdPara2.push_back(prodP2);
    ProdPara3.push_back(prodP3);
    EnergyMean.push_back(EnergyM);
    EnergySigma.push_back(EnergyS);
    Survival.push_back(SurvivalRate);
  }

  void Print(double all = -1) const 
  {
    std::cout<<" name BT:"<<BT<<" F:"<<F<<" "<< std::setprecision(8) << BT_ids[0]*100 + BT_ids[1] + 0.01*BT_ids[2] + 0.0001*BT_ids[3] << std::endl;
    std::cout<<" '-> N#"<<TargetCm.size()<<std::endl;
    for(unsigned int i=0;i<TargetCm.size();++i)
      {
	if(all < 0 || TMath::Abs(TargetCm[i]-all)<1.e-3)
	  {
	    std::cout<<" '---> Tcm:"<<TargetCm[i]<<" gCm:"<<TargetGramCm[i]<<std::endl;
	    std::cout<<" '---> Energy:"<<EnergyMean[i]<<" +/- "<<EnergySigma[i]<<" "<<Survival[i]<<std::endl;
	    std::cout<<" '---> Int:"<<ProdFrag[i]<<" C1:"<<ProdPara1[i]<<" C2:"<<ProdPara2[i]<<" C3:"<<ProdPara3[i]<<std::endl;
	  }
      }
  }
};

struct dataMax
{
  int idI;
  int idJ;
  double alpha;
  double beta;

  double beam;
  double target;
  double fragment;
  int cmtarget;
  double intensity;
  double contamination;
  double energyMean;
  double energySigma;
  double survivalRate;
  double crossSectionHyp;
};


class CompdataMax
{
public:
  bool operator() (const dataMax& a, const dataMax& b) const
  {
    if(a.idI == b.idI && a.idJ == b.idJ)
      {
	if(a.beam < b.beam)
	  return true;
	else
	  {
	    if(a.beam == b.beam)
	      {
		if(a.target < b.target)
		  return true;
		else
		  {
		    if(a.target == b.target)
		      {
			if(a.fragment < b.fragment)
			  return true;
			else
			  {
			    if(a.fragment == b.fragment)
			      {
				return a.cmtarget < b.cmtarget;
			      }
			    else
			      return false;
			  }
		      }
		    else
		      return false;
		  }
	      }
	    else
	      return false;
	  }
      }
    else
      {
	if(a.idI == b.idI)
	  return a.idJ < b.idJ;
	else
	  return a.idI < b.idI;
      }
      
  }
};


class nameBTtuple
{
public:
  std::string nameB;
  std::string nameT;
  double cross_section;

  // nameBTtuple() = default;
  // nameBTtuple(const nameBTtuple&) = default;
  // nameBTtuple(nameBTtuple&&) = default;
  // nameBTtuple& operator=(nameBTtuple&&) = default;
  
  //  nameBTtuple& operator=(const nameBTtuple& other ) 
  //  {
  //    nameBTtuple temp(other);
  //    std::swap(*this, temp); 
  //    return *this;
  //  }

};

class CompBTTuple
{
public :
  bool operator() (const nameBTtuple& a, const nameBTtuple& b) const
  {
    if(a.nameB!=b.nameB)
      return a.nameB < b.nameB;
    else
      return a.nameT < b.nameT;
  }
};



#endif // FRSANA_H
