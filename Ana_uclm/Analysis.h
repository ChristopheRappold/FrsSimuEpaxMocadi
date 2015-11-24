#ifndef FRSANA_H
#define FRSANA_H

#include <vector>
#include <list>
#include <string>
#include <iostream>

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

  void Print() const 
  {
    std::cout<<" name BT:"<<BT<<" F:"<<F<<" "<< std::setprecision(8) << BT_ids[0]*100 + BT_ids[1] + 0.01*BT_ids[2] + 0.0001*BT_ids[3] << std::endl;
    std::cout<<" '-> N#"<<TargetCm.size()<<std::endl;
    for(unsigned int i=0;i<TargetCm.size();++i)
      {
	std::cout<<" '---> Tcm:"<<TargetCm[i]<<" gCm:"<<TargetGramCm[i]<<std::endl;
	std::cout<<" '---> Energy:"<<EnergyMean[i]<<" +/- "<<EnergySigma[i]<<" "<<Survival[i]<<std::endl;
	std::cout<<" '---> Int:"<<ProdFrag[i]<<" C1:"<<ProdPara1[i]<<" C2:"<<ProdPara2[i]<<" C3:"<<ProdPara3[i]<<std::endl;
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



struct nameBTtuple
{
  std::string nameB;
  std::string nameT;
  double cross_section;
};

class CompBTTuple
{
public :
  bool operator() (const nameBTtuple& a, const nameBTtuple& b) const
  {
    return a.nameB < b.nameB;
  }
};



#endif // FRSANA_H
