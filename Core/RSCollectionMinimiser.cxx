#include"RSCollectionMinimiser.h"
#include<iostream>

// ROOT
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"

//////////////////////////////////////////////////////////////////////////////// 
RS::CollectionMinimiser::CollectionMinimiser(){
}

////////////////////////////////////////////////////////////////////////////////
RS::CollectionMinimiser::~CollectionMinimiser(){
}

////////////////////////////////////////////////////////////////////////////////
double RS::CollectionMinimiser::Minimise(std::vector<TGraphErrors*> Collection, std::vector< std::vector<TGraph*> > Bases, std::vector<double>& param, std::vector<double>& err, double lmin, double lmax){
  m_Collection = Collection;
  m_Bases = Bases;
  if(lmin<lmax){
    m_min= lmin;
    m_max= lmax;
  }

  else{
    m_min= lmax;
    m_max= lmin;
  }

  // Sanity check of the Bases:
  size_t mysize = m_Collection.size();
  if(mysize!= m_Bases.size()){
    std::cout << "ERROR: The collection size and the bases size provided to the Collection Minimiser are different" << std::endl;
    exit(1);
  }


  size_t base_size = 0 ;
  for(size_t i = 0 ; i < mysize; i++){
    if(i == 0) 
      base_size = m_Bases[i].size();
    else{
      size_t current_base_size = m_Bases[i].size();
      if(current_base_size!=base_size){
        std::cout << "ERROR: Bases provided to the Collection minimiser are of different size" << std::endl;
        exit(1);
      }
      base_size = current_base_size;
    }
  }
  // Resize the parameter and error
  param.resize(base_size,1);
  err.resize(base_size,0);

  const char* minName ="Minuit2";const char* algoName="Combined";

  ROOT::Math::Minimizer* min =
    ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetValidError(true);
  min->SetMaxFunctionCalls(1e6) ;
  // create funciton wrapper for minimizer
  ROOT::Math::Functor f(this,&CollectionMinimiser::ComputeOverallChi2, base_size);
  min->SetFunction(f);
  min->SetPrecision(1e-20);
  double* parameter = new double[base_size];
  for(unsigned int i = 0 ; i < base_size ; i++){
    parameter[i] = param[i];
    char name[4];
    sprintf(name,"V%d",i+1);

      min->SetVariable(i,name,parameter[i],100);
 
  }
  // do the minimization
  min->Minimize();
  const double *xs = min->X();
  for(size_t i = 0 ; i < base_size ; i++){
    param[i] = xs[i];
    err[i] = min->Errors()[i];
  }

  m_Collection.clear();
  m_Bases.clear();



  return ComputeOverallChi2(xs);
}

////////////////////////////////////////////////////////////////////////////////
double RS::CollectionMinimiser::ComputeOverallChi2(const double* param){
  size_t mysize = m_Collection.size();
  double chi2 = 0;

  for(size_t i =0; i < mysize ; i++){
    chi2+=ComputeOneChi2(i,param);
  }

  return chi2;
}

////////////////////////////////////////////////////////////////////////////////
double RS::CollectionMinimiser::ComputeOneChi2(unsigned int i, const double* param){
  size_t base_size = m_Bases[i].size();
  double chi=0;
  double chi2=0;
  int n = m_Collection[i]->GetN();
  double* x =  m_Collection[i]->GetX();
  double* y = m_Collection[i]->GetY();
  double* err = m_Collection[i]->GetEY();
  
  for(int point = 0 ; point < n ; point++){
    if(x[point]>m_min && x[point]<m_max){
      double th = 0;
      
      for(size_t j = 0 ; j < base_size; j++){
        th+= param[j]*m_Bases[i][j]->Eval(x[point]);
      }
      
      if(err[point]>0)
        chi = (y[point]-th)/err[point];
      else
        chi=0;
      
      chi2+=chi*chi;
    }
  }

  return chi2;
}
