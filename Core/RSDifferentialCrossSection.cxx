// STL
#include<iostream>
#include<fstream>


// RS++
#include"RSDifferentialCrossSection.h"
using namespace RS;
////////////////////////////////////////////////////////////////////////////////
DifferentialCrossSection::DifferentialCrossSection(){
  m_g = NULL;
  m_modified = false;
}
////////////////////////////////////////////////////////////////////////////////
DifferentialCrossSection::~DifferentialCrossSection(){
  if(m_g)
    delete m_g;
}
////////////////////////////////////////////////////////////////////////////////
void DifferentialCrossSection::LoadFromASCII(std::string filename){
  std::ifstream file(filename.c_str());
  double x,y;
  m_x.clear() ; m_y.clear();
  while(file >> x >> y){
    m_x.push_back(x);
    m_y.push_back(y);
  }
}
////////////////////////////////////////////////////////////////////////////////
double DifferentialCrossSection::Integrate(double low, double up){
 
 if(low>up){
   double temp = low;
   low = up;
   up = temp;
 }

 double integral = 0 ;
 unsigned int mysize = m_y.size();
 double previous_x = 0;
 for(unsigned int i = 0 ; i < mysize ; i++){
   if(m_x[i]>low && m_x[i]< up){
    integral += m_y[i]*sin(m_x[i]*M_PI/180.)*2*M_PI*(m_x[i]-previous_x)*M_PI/180.;
   }
   previous_x = m_x[i];
 }
 
 return integral;
}
////////////////////////////////////////////////////////////////////////////////
void DifferentialCrossSection::Scale(const double& scale){
 unsigned int mysize = m_y.size();
 for(unsigned int i = 0 ; i < mysize ; i++){
   m_y[i]*=scale;
 }
  m_modified = true;
}
////////////////////////////////////////////////////////////////////////////////
double DifferentialCrossSection::Eval(const double& x) {
  return GetTGraph()->Eval(x);
}

////////////////////////////////////////////////////////////////////////////////
DifferentialCrossSection& DifferentialCrossSection::operator+=(DifferentialCrossSection& right){
 unsigned int mysize = m_y.size();
 for(unsigned int i = 0 ; i < mysize ; i++){
   m_y[i]+=right.Eval(m_x[i]);
 }
 m_modified = true;
 return *this;
}
////////////////////////////////////////////////////////////////////////////////
DifferentialCrossSection& DifferentialCrossSection::operator-=(DifferentialCrossSection& right){
 unsigned int mysize = m_y.size();
 for(unsigned int i = 0 ; i < mysize ; i++){
   m_y[i]-=right.Eval(m_x[i]);
 }
 m_modified = true;
 return *this;
}


////////////////////////////////////////////////////////////////////////////////
DifferentialCrossSection operator+(DifferentialCrossSection& left, DifferentialCrossSection& right){
  left+=right;
  return left;
}
////////////////////////////////////////////////////////////////////////////////
DifferentialCrossSection operator-(DifferentialCrossSection& left, DifferentialCrossSection& right){
  left-=right;
  return left;
}
////////////////////////////////////////////////////////////////////////////////
TGraph* DifferentialCrossSection::GetTGraph(){
  if(!m_g)
    m_g = new TGraph(m_x.size(),&m_x[0],&m_y[0]);
  
  else if(m_modified){
    delete m_g;
    m_g = new TGraph(m_x.size(),&m_x[0],&m_y[0]);
    m_modified = false;
  }

  return m_g;
}

