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


template <typename T>
vector<size_t> sort_indexes(const vector<T> &v)
{
  
  // initialize original index locations
  vector<size_t> idx(v.size());
  for (size_t i = 0; i != idx.size(); ++i) idx[i] = i;
  
  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),[&v](size_t i1, size_t i2) {return v[i1] < v[i2];});
  return idx;
}

class Compute
{
  double f_logistic(double x, double x0, double k) { return 1./(1.+ TMath::Exp(-k*(x-x0)));};
public:

  double ComputeWeigth(double alpha, double beta, int Ntype)
  {
    switch(Ntype)
      {
      case 0 :
	return TMath::Sqrt(3./4. - alpha*alpha - beta*beta) ;
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
	return (alpha*alpha + beta*beta <= 3./4.);
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
	else if(type == 3)
	  {
	    In[i] = f_logistic(In[i],0.5*(Range[i][1]+Range[i][0]),-TMath::Log(1./(0.73*(Range[i][1]-0.5*(Range[i][1]+Range[i][0])))-1));
	  }
	else if(type == 4)
	  {
	    In[i] = f_logistic(In[i],Range[i][0],1./(0.5*(Range[i][2]-Range[i][1])));
	  }
	if(h_hist!=nullptr)
	  h_hist->Fill(6+i,In[i]);
      }
  }
  
  //TGeoManager *geom; 
  TGeoElementTable *tableRN; 
  std::map<std::string,double> cache_NucleiMass;

  Compute():tableRN(nullptr)//:geom(new TGeoManager("geom","radionuclides")),tableRN(geom->GetElementTable())
  { }
  
  double Brho(int Z, int A, double Ek)
  {
    if(gGeoManager==nullptr)
      {
	gGeoManager = new TGeoManager("geom","radionuclides");
      }
    if(tableRN==nullptr)
      tableRN = gGeoManager->GetElementTable();
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
    //std::cout<<A<<" "<<Z<<" "<<Mass<<" "<<p<<std::endl;
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

class HypDataProd
{
public:
  std::string nameB;
  std::string nameT;
  double cross_section;
  double SecondTargetD;
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

class CompHypDataProd
{
public :
  bool operator() (const HypDataProd& a, const HypDataProd& b) const
  {
    if(a.nameB!=b.nameB)
      return a.nameB < b.nameB;
    else
      return a.nameT < b.nameT;
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
  std::vector<double> Trans;

  tuple_graph()
  {}
  tuple_graph(const std::string& Name, const std::string& NameF, const std::vector<int>& ids, double TCm, double TgCM, double prodF, double prodP1, double prodP2, double prodP3,const std::vector<double>& stage, double EnergyM = 0, double EnergyS = 0, double SurvivalRate = 0,double Transmition = 0):
    BT(Name),F(NameF),BT_ids(ids),TargetCm(1,TCm),TargetGramCm(1,TgCM),ProdFrag(1,prodF),ProdPara1(1,prodP1),ProdPara2(1,prodP2),ProdPara3(1,prodP3),ProdStage(1,stage),EnergyMean(1,EnergyM),EnergySigma(1,EnergyS),Survival(1,SurvivalRate),Trans(1,Transmition)
  {}
  tuple_graph(const std::string& Name,const std::string& NameF , const std::vector<int>& ids, double TCm, double TgCM, double prodF, double prodP1, double prodP2, double prodP3, double EnergyM = 0, double EnergyS = 0, double SurvivalRate = 0,double Transmition = 0):
    BT(Name),F(NameF),BT_ids(ids),TargetCm(1,TCm),TargetGramCm(1,TgCM),ProdFrag(1,prodF),ProdPara1(1,prodP1),ProdPara2(1,prodP2),ProdPara3(1,prodP3),EnergyMean(1,EnergyM),EnergySigma(1,EnergyS),Survival(1,SurvivalRate),Trans(1,Transmition)
  {}
  ~tuple_graph() {}
  tuple_graph(const tuple_graph& t):BT(t.BT),F(t.F),BT_ids(t.BT_ids),TargetCm(t.TargetCm),TargetGramCm(t.TargetGramCm),ProdFrag(t.ProdFrag),ProdPara1(t.ProdPara1),ProdPara2(t.ProdPara1),ProdPara3(t.ProdPara3),ProdStage(t.ProdStage),EnergyMean(t.EnergyMean),EnergySigma(t.EnergySigma),Survival(t.Survival),Trans(t.Trans)
  {}
  void AddValue(double TCm, double TgCM, double prodF, double prodP1, double prodP2, double prodP3,const std::vector<double>& stage, double EnergyM = 0, double EnergyS = 0, double SurvivalRate = 0, double Transmition = 0)
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
    Trans.push_back(Transmition);
  }
  void AddValue(double TCm, double TgCM, double prodF, double prodP1, double prodP2, double prodP3, double EnergyM = 0, double EnergyS = 0, double SurvivalRate = 0, double Transmition=0)
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
    Trans.push_back(Transmition);
  }

  void Print(double all = -1) const 
  {
    std::cout<<" name BT:"<<BT<<" F:"<<F<<" "<< std::setprecision(8) << BT_ids[0]*100 + BT_ids[1] + 0.01*BT_ids[2] + 0.0001*BT_ids[3] << std::endl;
    std::cout<<" '-> N#"<<TargetCm.size()<<std::endl;

    std::vector<size_t> sorting = sort_indexes(TargetCm);

    for(unsigned int id=0;id<sorting.size();++id)
      {
	int i = sorting[id];
	if(all < 0 || TMath::Abs(TargetCm[i]-all)<1.e-3)
	  {
	    std::cout<<" '---> Tcm:"<<TargetCm[i]<<" gCm:"<<TargetGramCm[i]<<std::endl;
	    std::cout<<" '---> Energy:"<<EnergyMean[i]<<" +/- "<<EnergySigma[i]<<" "<<Survival[i]<<std::endl;
	    std::cout<<" '---> Int:"<<ProdFrag[i]<<" C1:"<<ProdPara1[i]<<" C2:"<<ProdPara2[i]<<" C3:"<<ProdPara3[i]<<std::endl;
	  }
      }
  }
};

struct Interpole{

  static constexpr double size_interval = 0.5;
  
  auto Inter (const tuple_graph& t)
  {
    tuple_graph t1(t);
    std::vector<size_t> idx_sorted = sort_indexes(t1.TargetCm);
    auto f = [] (const std::vector<double>& v, const std::vector<size_t> order) -> auto {
      std::vector<std::tuple<size_t,size_t> > res; 
      for(auto it1 = order.cbegin(), it2 = it1+1, it_end = order.cend(); it2 != it_end ;++it1, ++it2)
	{
	  if (TMath::Abs(v[*it1]-v[*it2])>size_interval)
	    {
	      res.push_back(std::make_tuple(*it1,*it2));
	    }
	}
      return res;
    };

    std::vector<std::tuple<size_t,size_t> > Intervals = f(t1.TargetCm,idx_sorted);
    // std::cout<<" Intervals :"<<Intervals.size();
    // for(const auto& interval : Intervals)
    //   std::cout<<" --- : "<<std::get<0>(interval)<<" "<<std::get<1>(interval)<<std::endl;
    
    for( const auto& interval : Intervals)
      {
	size_t id1 = std::get<0>(interval), id2 = std::get<1>(interval);
	double dist = t1.TargetCm[id2]-t1.TargetCm[id1];
	//std::cout<<" --- :"<<id1<<" "<<id2<<" "<<dist;

	if(dist<0)
	  std::swap(id1,id2);
	int nb_more = TMath::CeilNint(dist/size_interval)-1;
	//std::cout<<" "<<nb_more<<std::endl;
	
	std::vector<double> Alldelta;
	Alldelta.push_back((t1.TargetGramCm[id2]/t1.TargetCm[id2]   - t1.TargetGramCm[id1]/t1.TargetCm[id1])  /dist);
	Alldelta.push_back((t1.ProdFrag[id2]/(1.-t1.Survival[id2])  - t1.ProdFrag[id1]/(1.-t1.Survival[id1])) /dist);
	Alldelta.push_back((t1.ProdPara1[id2]/(1.-t1.Survival[id2]) - t1.ProdPara1[id1]/(1.-t1.Survival[id1]))/dist);
	Alldelta.push_back((t1.ProdPara2[id2]/(1.-t1.Survival[id2]) - t1.ProdPara2[id1]/(1.-t1.Survival[id1]))/dist);
	Alldelta.push_back((t1.ProdPara3[id2]/(1.-t1.Survival[id2]) - t1.ProdPara3[id1]/(1.-t1.Survival[id1]))/dist);
	Alldelta.push_back((t1.EnergyMean[id2]			    - t1.EnergyMean[id1])		      /dist);
	Alldelta.push_back((t1.EnergySigma[id2]			    - t1.EnergySigma[id1])		      /dist);
	Alldelta.push_back((t1.Survival[id2]			    - t1.Survival[id1])			      /dist);

	std::vector<double> Allini_b;
	Allini_b.push_back(t1.TargetGramCm[id1]/t1.TargetCm[id1]);  
	Allini_b.push_back(t1.ProdFrag[id1]/(1.-t1.Survival[id1])); 
	Allini_b.push_back(t1.ProdPara1[id1]/(1.-t1.Survival[id1]));
	Allini_b.push_back(t1.ProdPara2[id1]/(1.-t1.Survival[id1]));
	Allini_b.push_back(t1.ProdPara3[id1]/(1.-t1.Survival[id1]));
	Allini_b.push_back(t1.EnergyMean[id1]);			    
	Allini_b.push_back(t1.EnergySigma[id1]);		    
	Allini_b.push_back(t1.Survival[id1]);			  

	double ini_a = t1.TargetCm[id1];

	std::vector<double> res(Alldelta.size());
	auto Dointer = [&Alldelta,&Allini_b,&ini_a,&res] (int index, double x) {
	  res[index] = Allini_b[index] + Alldelta[index]*(x-ini_a);
	  //std::cout<<" Int: ["<<index<<"] "<<x<<" : "<<res[index]<<" | "<<Alldelta[index]<<" "<<Allini_b[index]<<" "<<ini_a<<std::endl;
	};
	
	for(int i = 0; i < nb_more; ++i)
	  {
	    double temp_Cm = ini_a+size_interval*(i+1.);
	    for(size_t j = 0; j<res.size();++j)
	      Dointer(j,temp_Cm);

	    t1.AddValue(temp_Cm,res[0]*temp_Cm,res[1]*(1.-res[7]),res[2]*(1.-res[7]),res[3]*(1.-res[7]),res[4]*(1.-res[7]),res[5],res[6],res[7]);
	  }	
      }

    return t1;
  };
  
};

struct Histo
{
  bool Dohisto;
  TH2F* h_hist;
  TH2F* h_histNameBeamPrimary;
  TH2F* h_histNameBeamSecondary;
  TH2F* h_histCrossSection;
  TH2F* h_histEnergy;
  TH2F* h_histCmTarget;
  TH2F* h_histNameTarget;
  TH2F* h_Internal;
};

class Regularize
{
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


public:
  double Tminmax[2];
  double Iminmax[2];
  double Cminmax[2];
  double ICminmax[2];
  double Energyminmax [2];
  double Survivalminmax [2];
  double Transminmax [2];
  double CXminmax[2];
  
  double Tmeansigma [2];
  double Imeansigma [2];
  double Cmeansigma [2];
  double ICmeansigma [2];
  double Energymeansigma [2];
  double Survivalmeansigma [2];
  double Transmeansigma [2];
  double CXmeansigma [2];
  
  double Tquantile [3];
  double Iquantile [3];
  double Cquantile [3];
  double ICquantile [3];
  double Energyquantile [3];
  double Survivalquantile [3];
  double Transquantile [3];
  double CXquantile [3];
  std::vector<std::vector<double*> > RangeAll = {{Tminmax,Iminmax,Cminmax,Energyminmax,CXminmax,Survivalminmax,Transminmax},
						 {Tmeansigma,Imeansigma,Cmeansigma,Energymeansigma,CXmeansigma,Survivalmeansigma,Transmeansigma},
						 {Tquantile,Iquantile,Cquantile,Energyquantile,CXquantile,Survivalquantile,Transquantile},
						 {Tminmax,Iminmax,Cminmax,Energyminmax,CXminmax,Survivalminmax,Transminmax},
						 {Tmeansigma,Imeansigma,Cmeansigma,Energymeansigma,CXmeansigma,Survivalmeansigma,Transmeansigma},
						 {Tquantile,Iquantile,Cquantile,Energyquantile,CXquantile,Survivalquantile,Transquantile}  };
  
  Regularize(std::vector<double>& targetLength,std::vector<double>& intensity, std::vector<double>&contamination, std::vector<double>& IntCon, std::vector<double>& EnergyM, std::vector<double>& Survival, std::vector<double>& Trans,std::vector<double>& CX, const std::string& ICtype, double IntensityMax)
  {

    Tminmax[0] = *std::min_element(targetLength.begin(),targetLength.end());
    Iminmax[0] = *std::min_element(intensity.begin(),intensity.end());
    Cminmax[0] = *std::min_element(contamination.begin(),contamination.end());
    ICminmax[0] = *std::min_element(IntCon.begin(),IntCon.end());
    Energyminmax[0]  =*std::min_element(EnergyM.begin(),EnergyM.end());
    Survivalminmax[0] =*std::min_element(Survival.begin(),Survival.end());
    Transminmax[0] =*std::min_element(Trans.begin(),Trans.end());
    CXminmax[0] =*std::min_element(CX.begin(),CX.end());

    
    Tminmax[1] = *std::max_element(targetLength.begin(),targetLength.end());
    Iminmax[1] = *std::max_element(intensity.begin(),intensity.end());
    Cminmax[1] = *std::max_element(contamination.begin(),contamination.end());
    ICminmax[1] = *std::max_element(IntCon.begin(),IntCon.end());
    Energyminmax[1]  =*std::max_element(EnergyM.begin(),EnergyM.end());
    Survivalminmax[1] =*std::max_element(Survival.begin(),Survival.end());
    Transminmax[1] =*std::max_element(Trans.begin(),Trans.end());
    CXminmax[1] =*std::max_element(CX.begin(),CX.end());

    Iminmax[1] = Iminmax[1] > IntensityMax ? IntensityMax : Iminmax[1];

    computeMeanSigma(intensity,Imeansigma);
    computeMeanSigma(targetLength,Tmeansigma);
    computeMeanSigma(contamination,Cmeansigma);
    computeMeanSigma(IntCon,ICmeansigma);
    computeMeanSigma(EnergyM, Energymeansigma);
    computeMeanSigma(Survival, Survivalmeansigma);
    computeMeanSigma(Trans, Transmeansigma);
    computeMeanSigma(CX, CXmeansigma);

  
    computeQuantile(intensity,Iquantile);
    computeQuantile(targetLength,Tquantile);
    computeQuantile(contamination,Cquantile);
    computeQuantile(IntCon,ICquantile);
    computeQuantile(EnergyM, Energyquantile);
    computeQuantile(Survival, Survivalquantile);
    computeQuantile(Trans, Transquantile);
    computeQuantile(CX, CXquantile);

    if(ICtype == "I-C")
    {
      RangeAll[0][2] = ICminmax;
      RangeAll[1][2] = ICmeansigma;
      RangeAll[2][2] = ICquantile;
    }

  }
  


};

class RejectedBeamStruct
{
public:
  std::set<std::string> RejectedBeam;
  RejectedBeamStruct()
  {
    const std::map<int,int> StableBeam2 = {{1,1},{1,2},{2,3},{2,4},{3,6},{3,7},{4,9},{5,10},{5,11},{6,12},{6,13},{7,14},{7,15},{8,16},{8,17},{8,18},{9,19},{10,20},{10,21},{10,22},{11,23},{12,24},{12,25},{12,26},{13,27},{14,28},{14,29},{14,30},{15,31},{16,32},{16,33},{16,34},{16,36},{17,35},{17,37},{18,36},{18,38},{18,40}};

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
  }
  
};


class Minimizer
{
  Compute& compute;
  //const Regularize& reg;
  const std::string& ICtype;
  int type;
  const std::string& Mtype;
  int NtypeComp;
  const LambdaCrossSectionFunction StrangenessProdNorm;
  double minCX;
public :
  Minimizer(Compute& c,const std::string& _ICtype, int _type,const std::string& _Mtype, int _NtypeComp, double LambdaEnergy, double _minCX=0.5):compute(c),ICtype(_ICtype),type(_type),Mtype(_Mtype),NtypeComp(_NtypeComp),StrangenessProdNorm(LambdaCrossSectionFunction(LambdaEnergy)),minCX(_minCX)
  {
  }
  
  auto DoMin(double alpha, double beta, size_t idCmTarget, const tuple_graph& t1, const HypDataProd& HypData, Histo& histo, const Regularize& reg)
  {
    std::vector<int> tempIds(t1.BT_ids);
    
                             
    //std::cout<<" -primary:"<<t1.F<<" "<<it_FBT->first<<std::endl;
    //std::cout<<" idCm :"<<idCmTarget<<std::endl;
    double tempT = t1.TargetGramCm[idCmTarget];
    double tempI = t1.ProdFrag[idCmTarget];
    if(tempI>reg.Iminmax[1])
      tempI=reg.Iminmax[1];
    double tempEpar[2] = {t1.EnergyMean[idCmTarget],t1.EnergySigma[idCmTarget]}; 
    
    double tempICX = tempI*4.*HypData.SecondTargetD*0.2*HypData.cross_section*6.02214129e-07*StrangenessProdNorm(tempEpar[0]); // 6.02409638554217e-07 = 1e-30/1.66e-24
    double tempC = t1.ProdPara1[idCmTarget]+t1.ProdPara2[idCmTarget]+t1.ProdPara3[idCmTarget];
    //double tempE = (2.e3-tempEpar[0])*(2.e3-tempEpar[0])/tempEpar[1]/tempEpar[1];
    if(ICtype == "I-C")
      tempC = tempI - tempC;
    // if(t1.BT_ids[0] == 14 && t1.BT_ids[1] == 7 && t1.BT_ids[2] == 9 && t1.BT_ids[3] == 4 && t1.BT_ids[4] == 12 && t1.BT_ids[5] == 7 && TMath::Abs(t1.TargetCm[idCmTarget]-5.)<1.e-3 && it_secondaryReaction->nameT=="C12")
    //   {
    //     std::cout<<"Test: "<<it_FBT->first<<std::endl;
    //     t1.Print(5.);
    //     std::cout<<" Secondary: "<<it_secondaryReaction->nameB<<" "<<it_secondaryReaction->nameT<<" "<<tempICX<<std::endl;
    //   }
    
    //h_Internal->Fill(5.,tempICX/640.);
    
    if(tempICX<minCX)
      return std::make_tuple(-2,0.,0.,tempIds,0.,0.,0.,0.,0.,0.);
    
    double Brho = compute.Brho(t1.BT_ids[5],t1.BT_ids[4],tempEpar[0]); 
    if(histo.Dohisto==false)
      histo.h_Internal->Fill(14.,Brho/10.);
    
    if(Brho>20.5)
      return std::make_tuple(-3,0.,0.,tempIds,0.,0.,0.,0.,0.,0.);
    //double tempE = 10*TMath::Gaus(2000,tempEpar[0],tempEpar[1]*5);
    
    //h_Internal->Fill(11.,(tempT-Tmeansigma[0])/(Tmeansigma[1]));
    //h_Internal->Fill(12.,(tempI-Imeansigma[0])/(Imeansigma[1]));
    //h_Internal->Fill(13.,(tempEpar[0]-Energymeansigma[0])/(Energymeansigma[1]));
    
    
    std::vector<double> Input = {tempT,tempI,tempC,tempEpar[0],tempICX};
    //std::vector<double> Input2 = {0,tempICX};
    
    compute.Invariant(Input, reg.RangeAll[type], type, histo.h_Internal);
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
      return std::make_tuple(-100,0.,0.,tempIds,0.,0.,0.,0.,0.,0.);
    if(res==-2)
      return std::make_tuple(-4,0.,0.,tempIds,0.,0.,0.,0.,0.,0.);
    //double tempG =  - 1./0.5808*1./0.54*1/1.10*outPar[0]*Input[0] + 1./320*1./0.32*outPar[1]*tempICX + outPar[2]*Input[3] ;//+ outPar[2]*Input[2] ;

    //double tempG =  - 1./0.5808*1./0.54*1/1.10*1./0.92*outPar[0]*Input[0] + (1./320*1./0.32*1./0.14*1./1.02*1./0.98)*outPar[1]*tempICX + 1./0.98*outPar[2]*((Input[3]-0.46)/(1.-0.46)) + 1./1.02*0.5*Input[1];//+ outPar[2]*Input[2]

    double tempG =  -outPar[0]*Input[0] + outPar[1]*Input[4] + outPar[2]*Input[3] + 0.5*Input[1];//+ outPar[2]*Input[2]


    //double tempG =  - 1./0.5808*outPar[0]*Input[0] + 1./0.56448*outPar[1]*Input[1] + outPar[2]*Input[3] + 1./5.4*Input2[1];
		  
    // h_Internal->Fill(0.,tempT);
    // h_Internal->Fill(1.,tempI);
    // h_Internal->Fill(2.,tempC);
    // h_Internal->Fill(3.,tempEpar[0]);
    // h_Internal->Fill(4.,tempG);
    if(histo.Dohisto==false)
      {
	histo.h_Internal->Fill(0.,Input[0]);
	histo.h_Internal->Fill(1.,Input[1]);
	histo.h_Internal->Fill(2.,Input[4]);
	histo.h_Internal->Fill(3.,Input[3]);
      }
    return std::make_tuple(0,tempG,t1.TargetCm[idCmTarget],tempIds,t1.ProdFrag[idCmTarget], t1.ProdPara1[idCmTarget]+t1.ProdPara2[idCmTarget]+t1.ProdPara3[idCmTarget],t1.EnergyMean[idCmTarget],t1.EnergySigma[idCmTarget], t1.Survival[idCmTarget],tempICX);
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


  auto CheckConditions(double alpha, double beta, size_t idCmTarget, const tuple_graph& t1, const HypDataProd& HypData, double IntensityMax)
  {
    //std::cout<<" -primary:"<<t1.F<<" "<<it_FBT->first<<std::endl;
    //std::cout<<" idCm :"<<idCmTarget<<std::endl;
    double tempT = t1.TargetGramCm[idCmTarget];
    double tempI = t1.ProdFrag[idCmTarget];
    if(tempI>IntensityMax)
      tempI=IntensityMax;
    double tempEpar[2] = {t1.EnergyMean[idCmTarget],t1.EnergySigma[idCmTarget]}; 
    
    double tempICX = tempI*4.*HypData.SecondTargetD*0.2*HypData.cross_section*6.02214129e-07*StrangenessProdNorm(tempEpar[0]); // 6.02409638554217e-07 = 1e-30/1.66e-24
    double tempC = t1.ProdPara1[idCmTarget]+t1.ProdPara2[idCmTarget]+t1.ProdPara3[idCmTarget];

    if(ICtype == "I-C")
      tempC = tempI - tempC;
    // if(t1.BT_ids[0] == 14 && t1.BT_ids[1] == 7 && t1.BT_ids[2] == 9 && t1.BT_ids[3] == 4 && t1.BT_ids[4] == 12 && t1.BT_ids[5] == 7 && TMath::Abs(t1.TargetCm[idCmTarget]-5.)<1.e-3 && it_secondaryReaction->nameT=="C12")
    //   {
    //     std::cout<<"Test: "<<it_FBT->first<<std::endl;
    //     t1.Print(5.);
    //     std::cout<<" Secondary: "<<it_secondaryReaction->nameB<<" "<<it_secondaryReaction->nameT<<" "<<tempICX<<std::endl;
    //   }
    
    if(tempICX<minCX)
      return std::make_tuple(-2,0.,tempICX);
    
    double Brho = compute.Brho(t1.BT_ids[5],t1.BT_ids[4],tempEpar[0]); 
    
    if(Brho>20.5)
      return std::make_tuple(-3,0.,tempICX);

    std::vector<double> Input = {tempT,tempI,tempC,tempEpar[0]};

    std::vector<double> inPar(2);
    inPar[0] = alpha;
    inPar[1] = beta;

    std::vector<double> outPar(3);
    int res = compute.setPermutation(Mtype,inPar,outPar,NtypeComp);
    if(res==-1)
      return std::make_tuple(-100,0.,tempICX);
    if(res==-2)
      return std::make_tuple(-4,0.,tempICX);
    //double tempG =  - 1./0.5808*1./0.54*1/1.10*outPar[0]*Input[0] + 1./320*1./0.32*outPar[1]*tempICX + outPar[2]*Input[3] ;//+ outPar[2]*Input[2] ;
    double tempG =  outPar[0]*Input[0] + outPar[1]*tempICX + outPar[2]*Input[3] + 0.5*Input[1];//+ outPar[2]*Input[2]

    return std::make_tuple(0,tempG,tempICX);
	        
  }


};

struct dataMax
{
  size_t idI;
  size_t idJ;
  double alpha;
  double beta;

  double beam;
  double target;
  double fragment;
  double cmtarget;
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





#endif // FRSANA_H
