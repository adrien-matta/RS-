// STL
#include<iostream>
#include<fstream>
#include<sstream>

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
}
////////////////////////////////////////////////////////////////////////////////
void DifferentialCrossSection::LoadFromASCII(std::string filename){
  std::ifstream file(filename.c_str());
  double x,y;
  double dump;
  m_x.clear() ; m_y.clear();
  while(file >> x >> y ){
    m_x.push_back(x);
    m_y.push_back(y);
  }
  m_modified = true;
}
////////////////////////////////////////////////////////////////////////////////
void DifferentialCrossSection::LoadFromTWOFNR(std::string filename){
  std::ifstream file(filename.c_str());

  // CM and LAB file have different formatting so find our which one:
  std::string LineBuffer;
  getline(file,LineBuffer);
  std::istringstream line(LineBuffer);
  double buffer;
  unsigned int mycount = 0 ;
  while(line >> buffer)
     mycount++;

  if(mycount==2){
    LoadFromASCII(filename);
    return;
  }
    
  file.close();
  file.open(filename.c_str());
  double x,y;
  double dump;
  m_x.clear() ; m_y.clear();
  while(file >> x >> y >> dump ){
    m_x.push_back(x);
    m_y.push_back(y);
  }

  m_modified = true;
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
  double previous_x = 0  ;
  if(mysize>1){
    if(m_x[1] < m_x[0])
      previous_x = 180;
  }
  for(unsigned int i = 0 ; i < mysize ; i++){
    if(m_x[i]>low && m_x[i]< up){
      integral += m_y[i]*sin(m_x[i]*M_PI/180.)*2*M_PI*fabs(m_x[i]-previous_x)*M_PI/180.;
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
  if(mysize==0) {
    mysize = right.m_x.size();
    m_x.resize(mysize,0);
    for(unsigned int i = 0 ; i < mysize ; i++){
      m_x[i]=right.m_x[i]; 
    }
    m_y.resize(mysize,0);
  }

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
DifferentialCrossSection RS::operator+(DifferentialCrossSection& left, DifferentialCrossSection& right){
  left+=right;
  return left;
}
////////////////////////////////////////////////////////////////////////////////
DifferentialCrossSection RS::operator-(DifferentialCrossSection& left, DifferentialCrossSection& right){
  left-=right;
  return left;
}
////////////////////////////////////////////////////////////////////////////////
TGraph* DifferentialCrossSection::GetTGraph(){
  if(!m_g)
    m_g = new TGraph(m_x.size(),&m_x[0],&m_y[0]);

  else if(m_modified){
    //delete m_g;
    m_g = new TGraph(m_x.size(),&m_x[0],&m_y[0]);
    m_modified = false;
  }

  return m_g;
}

