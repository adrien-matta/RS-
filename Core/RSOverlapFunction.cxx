#include"RSOverlapFunction.h"
#include<fstream>
#include<sstream>
#include<iostream>
#include<algorithm>

// ROOT
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TGraphErrors.h"
#include "TF1.h"
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
RS::OverlapFunction::OverlapFunction(){
  m_Graph=0;
  m_modified=false;
}
////////////////////////////////////////////////////////////////////////////////
RS::OverlapFunction::~OverlapFunction(){
}

////////////////////////////////////////////////////////////////////////////////
RS::OverlapFunction::OverlapFunction(std::string filename){
  LoadFromFile(filename);
}
////////////////////////////////////////////////////////////////////////////////
void RS::OverlapFunction::LoadFromFile(std::string filename){

  m_OverlapFunction_R.clear();
  m_OverlapFunction_OL.clear();
  m_OverlapFunction_E.clear();

  std::ifstream infile(filename.c_str());
  if(!infile.is_open()){
    std::cout << "Error canot open overlap file " << filename << std::endl;
    exit(1);
  }

  // Analyse the file to see if there is errors
  std::string Line,Buffer;
  std::getline(infile,Line);

  int column=0;
  std::istringstream iss(Line);
  while(iss>>Buffer)
    column++;

  double r,OL,e;
  if(column==2)
    while(infile >> r >> OL){

      m_OverlapFunction_R.push_back(r);
      m_OverlapFunction_OL.push_back(OL);
      m_OverlapFunction_E.push_back(0);
    }

  else if(column==3)
    while(infile >> r >> OL >> e){

      m_OverlapFunction_R.push_back(r);
      m_OverlapFunction_OL.push_back(OL);
      m_OverlapFunction_E.push_back(e);
    }

  else{
    std::cout << "Error : Wrong number of column in Overlap File : " << filename << std::endl;
    exit(1);
  }


  Sort();
}
////////////////////////////////////////////////////////////////////////////////
void RS::OverlapFunction::Sort(){
  std::vector<double> tempR = m_OverlapFunction_R;
  std::vector<double> tempR2 = m_OverlapFunction_R;

  std::vector<double> tempOL = m_OverlapFunction_OL;
  std::vector<double> tempE = m_OverlapFunction_E;

  sort(tempR.begin(),tempR.end());

  size_t mysize = tempR.size();
  for(size_t i = 0 ; i < mysize ; i++){
    for(size_t j = 0 ; j < mysize ; j++){
      if(tempR2[j] == tempR[i]){ 
        m_OverlapFunction_R[i] = tempR2[j];
        m_OverlapFunction_OL[i] = tempOL[j];
        m_OverlapFunction_E[i] = tempE[j];
      }
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
void RS::OverlapFunction::Normalise(double norm){
  double inte = Integrate();
  size_t mysize = m_OverlapFunction_OL.size();

  for(size_t i = 0 ; i < mysize ; i++){
    m_OverlapFunction_OL[i]=m_OverlapFunction_OL[i]*norm/inte;
  }
  m_modified=true;
}
////////////////////////////////////////////////////////////////////////////////
void RS::OverlapFunction::Scale(double factor){
  size_t mysize = m_OverlapFunction_OL.size();
  
  double efactor = factor;
  if(efactor<0)
    efactor=sqrt(factor*factor);

  for(size_t i = 0 ; i < mysize ; i++){
    m_OverlapFunction_OL[i]*=factor;
    m_OverlapFunction_E[i]*=efactor;
  }
  m_modified=true;
}
////////////////////////////////////////////////////////////////////////////////
void RS::OverlapFunction::Divide(TGraph* g){
  size_t mysize = m_OverlapFunction_OL.size();

  for(size_t i = 0 ; i < mysize ; i++){
    m_OverlapFunction_OL[i]=m_OverlapFunction_OL[i]/g->Eval(m_OverlapFunction_R[i]);
  }
  m_modified=true;
}
////////////////////////////////////////////////////////////////////////////////
void RS::OverlapFunction::Multiply(TF1* f){
  size_t mysize = m_OverlapFunction_OL.size();

  for(size_t i = 0 ; i < mysize ; i++){
    m_OverlapFunction_OL[i]*=f->Eval(m_OverlapFunction_R[i]);
  }
  m_modified=true;
  
}

////////////////////////////////////////////////////////////////////////////////
void RS::OverlapFunction::Add(TGraph* g,double f){
  size_t mysize = m_OverlapFunction_OL.size();

  for(size_t i = 0 ; i < mysize ; i++){
    m_OverlapFunction_OL[i]=m_OverlapFunction_OL[i]+f*g->Eval(m_OverlapFunction_R[i]);
  }
  m_modified =true;
}
////////////////////////////////////////////////////////////////////////////////
double RS::OverlapFunction::Integrate(double min, double max){
  if(min < 0)
    min = 0;


  double inte=0;
  double step=0;


  size_t mysize = m_OverlapFunction_OL.size();

  if(mysize > 1)
    step = m_OverlapFunction_R[1]-m_OverlapFunction_R[0];

  for(size_t i = 0 ; i < mysize ; i++){
    if(i<mysize-1)
      step = m_OverlapFunction_R[i]-m_OverlapFunction_R[i+1]; 

    inte += m_OverlapFunction_OL[i]*step;

  }

  return inte;
}
////////////////////////////////////////////////////////////////////////////////
TGraph* RS::OverlapFunction::GetTGraph(){

  if(m_Graph==NULL || m_modified){ 

    double sum = 0 ;

    size_t mysize = m_OverlapFunction_R.size();
    for(size_t i = 0 ; i < mysize ; i++){
      sum += m_OverlapFunction_E[i]; 
    }
    if(sum==0)
      m_Graph = new TGraph(mysize,&m_OverlapFunction_R[0],&m_OverlapFunction_OL[0]);

    else
      m_Graph = new TGraphErrors(mysize,&m_OverlapFunction_R[0],&m_OverlapFunction_OL[0],0,&m_OverlapFunction_E[0]);

    m_modified=false;
  }

  return m_Graph;

}
////////////////////////////////////////////////////////////////////////////////
double RS::OverlapFunction::Eval(double R){

  return GetTGraph()->Eval(R);
}
////////////////////////////////////////////////////////////////////////////////
RS::OverlapFunction& RS::OverlapFunction::operator+=(OverlapFunction& right){
  unsigned int mysize = m_OverlapFunction_OL.size();
  if(mysize==0) {
    mysize = right.m_OverlapFunction_R.size();
    m_OverlapFunction_R.resize(mysize,0);
    m_OverlapFunction_E.resize(mysize,0);
    for(unsigned int i = 0 ; i < mysize ; i++){
      m_OverlapFunction_R[i]=right.m_OverlapFunction_R[i];
      m_OverlapFunction_E[i]=right.m_OverlapFunction_E[i];
    }
    m_OverlapFunction_OL.resize(mysize,0);
  }

  for(unsigned int i = 0 ; i < mysize ; i++){
    m_OverlapFunction_OL[i]+=right.Eval(m_OverlapFunction_R[i]);
  }
  m_modified = true;
  return *this;
}
////////////////////////////////////////////////////////////////////////////////
RS::OverlapFunction& RS::OverlapFunction::operator-=(OverlapFunction& right){
  unsigned int mysize = m_OverlapFunction_OL.size();
  for(unsigned int i = 0 ; i < mysize ; i++){
    m_OverlapFunction_OL[i]-=right.Eval(m_OverlapFunction_R[i]);
  }
  m_modified = true;
  return *this;
}



////////////////////////////////////////////////////////////////////////////////
RS::OverlapFunction RS::OverlapFunction::FitWithOverlap(std::vector<RS::OverlapFunction> Th, std::vector<double>& param,std::vector<double>& err, double lmin,double lmax){
  // Make a local copy of the function
  if(lmin<lmax){
    m_min= lmin;
    m_max= lmax;
  }

  else{
    m_min= lmax;
    m_max= lmin;
  }

  m_ThOriginal.clear();
  m_ThOriginal = Th;

  unsigned int mysize = Th.size();
  param.resize(mysize,0);
  err.resize(mysize,0);

  const char* minName ="Minuit2";const char* algoName="Fumili2";

  ROOT::Math::Minimizer* min =
  ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetValidError(true);
  min->SetMaxFunctionCalls(1e6) ;
  // create funciton wrapper for minimizer
  // a IMultiGenFunction type
  ROOT::Math::Functor f(this,&OverlapFunction::ComputeChi2Fit, mysize);
  min->SetFunction(f);
  double* parameter = new double[mysize];
  for(unsigned int i = 0 ; i < mysize ; i++){
    parameter[i] = param[i];
    char name[4];
    sprintf(name,"V%d",i+1);
    min->SetLimitedVariable(i,name,parameter[i],0.1,-1000,1000);
  }
  // do the minimization
  min->Minimize();
  const double *xs = min->X();

  // Return the Fitted CS
  OverlapFunction OL;
  for(unsigned int i = 0 ; i < mysize ; i++){
    param[i] = xs[i];
    err[i] = min->Errors()[i]; 
    Th[i].Scale(xs[i]);
    OL+=Th[i];
  }

  // Clean up
  m_ThOriginal.clear();
  m_ThCurrent.clear(); 

  return OL;

}
////////////////////////////////////////////////////////////////////////////////
double RS::OverlapFunction::ComputeChi2(OverlapFunction& OF){
  TGraph* th = OF.GetTGraph();
  double chi2 = 0;
  double chi ;
  int mysize = m_OverlapFunction_R.size();
  for(int i = 0 ; i < mysize ; i++){
    if(m_OverlapFunction_E[i]!=0){
      if((m_min<0 || m_OverlapFunction_R[i]>m_min) && (m_max<0 || m_OverlapFunction_R[i]<m_max)){
        chi =  (m_OverlapFunction_OL[i]-th->Eval(m_OverlapFunction_R[i]))/m_OverlapFunction_E[i];
        chi2+=chi*chi;
      }
    }
  }
  if(chi2!=chi2)
    return 1e8;
  return chi2; 
}
////////////////////////////////////////////////////////////////////////////////
double RS::OverlapFunction::ComputeChi2Fit(const double* parameter){
  unsigned int mysize = m_ThOriginal.size();
  m_ThCurrent.resize(mysize);
  OverlapFunction OF;
  for(unsigned int i = 0 ; i < mysize ; i++){
    m_ThCurrent[i]= m_ThOriginal[i];
    m_ThCurrent[i].Scale(parameter[i]);
    OF+= m_ThCurrent[i];
  }
  return ComputeChi2(OF); 
}
