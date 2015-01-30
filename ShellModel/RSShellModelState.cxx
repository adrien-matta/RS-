#include"RSShellModelState.h"

#include<iostream>
#include<sstream>
#include<cmath>
using namespace RS;
////////////////////////////////////////////////////////////////////////////////
ShellModelState::ShellModelState(){
}

////////////////////////////////////////////////////////////////////////////////
ShellModelState::ShellModelState(double E, double J,int P,int order){
  m_E = E;
  m_J = J;
  m_P = P;
  m_order = order;
}

////////////////////////////////////////////////////////////////////////////////
ShellModelState::~ShellModelState(){
}

////////////////////////////////////////////////////////////////////////////////
double ShellModelState::GetEnergy(){
  return m_E;
}
////////////////////////////////////////////////////////////////////////////////
void ShellModelState::SetEnergy(double Energy){
  m_E = Energy;
}
////////////////////////////////////////////////////////////////////////////////
double ShellModelState::GetJ(){
  return m_J;
}

////////////////////////////////////////////////////////////////////////////////
int ShellModelState::GetParity(){
  return m_P;
}

////////////////////////////////////////////////////////////////////////////////
int ShellModelState::GetOrder(){
  return m_order;
}

////////////////////////////////////////////////////////////////////////////////
void ShellModelState::AddOrbital(unsigned int n, unsigned int l, double j, double s){
  m_n.push_back(n);
  m_l.push_back(l);
  m_j.push_back(j);
  m_s.push_back(s);
}

////////////////////////////////////////////////////////////////////////////////
unsigned int ShellModelState::GetNumberOfOrbital(){
  return m_n.size();
 }

////////////////////////////////////////////////////////////////////////////////
std::string ShellModelState::GetOrbitalString(unsigned int i){

  std::string shell="?";
  if(m_l[i]==0) shell ="s";
  else if(m_l[i]==1) shell ="p";
  else if(m_l[i]==2) shell ="d";
  else if(m_l[i]==3) shell ="f";
  else if(m_l[i]==4) shell ="g";
  else if(m_l[i]==5) shell ="h";
  else if(m_l[i]==6) shell ="i";
  else if(m_l[i]==7) shell ="j";
  else if(m_l[i]==8) shell ="k";

  std::ostringstream os;
  os << m_n[i] << shell << (m_j[i]*2) << "/2";

  return os.str();
}

////////////////////////////////////////////////////////////////////////////////
unsigned int ShellModelState::GetOrbitalN(unsigned int i){
  return m_n[i];
}

////////////////////////////////////////////////////////////////////////////////
unsigned int ShellModelState::GetOrbitalL(unsigned int i){
  return m_l[i];
}

////////////////////////////////////////////////////////////////////////////////
double ShellModelState::GetOrbitalJ(unsigned int i){
  return m_j[i];
}

////////////////////////////////////////////////////////////////////////////////
double ShellModelState::GetOrbitalS(unsigned int i){
  return m_s[i];
}

////////////////////////////////////////////////////////////////////////////////
unsigned int ShellModelState::GetMainOrbital(){
  unsigned int orbital = 1000;
  double s = -1;
  unsigned int mysize = m_n.size();
  
  for(unsigned int i = 0 ; i < mysize ; i++){
    if(m_s[i]>s){
       orbital = i;
       s = m_s[i];
    }
  }

  return orbital;

}

////////////////////////////////////////////////////////////////////////////////
double ShellModelState::GetSumOfSForL(unsigned int l){
  double sum = 0;
  unsigned int mysize = m_n.size();
  
  for(unsigned int i = 0 ; i < mysize ; i++){
    if(m_l[i]==l)
       sum+= m_s[i];
  }

  return sum;
}

////////////////////////////////////////////////////////////////////////////////
void ShellModelState::Print(){
  std::cout << "////// State E= " << m_E << " J= " << m_J << "(" << m_order << ")";
  if(m_P<0) std::cout << "-" << std::endl;
  else      std::cout << "+" << std::endl;

 unsigned int mysize = m_n.size();
 for(unsigned int i = 0 ; i < mysize ; i++){
   std::cout << GetOrbitalString(i) << " S= " << m_s[i] << std::endl;
 } 
}
////////////////////////////////////////////////////////////////////////////////
