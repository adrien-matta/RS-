#include"RSShellModelState.h"

#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<set>
#include <dirent.h>


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
  // Guarantee that the orbital list is by SF ordering importance (biggest SF first)

  if(!m_n.empty()){
    std::vector<unsigned int>::iterator it_n = m_n.begin();;
    std::vector<unsigned int>::iterator it_l = m_l.begin();;
    std::vector<double>::iterator it_j = m_j.begin();
    std::vector<double>::iterator it_s = m_s.begin();
    bool done = false;
    for(it_n = m_n.begin(); it_n!=m_n.end(); it_n++,it_l++,it_j++,it_s++){
      if(*it_s <= s){
        m_n.insert(it_n,n);
        m_l.insert(it_l,l);
        m_j.insert(it_j,j);
        m_s.insert(it_s,s);
        done = true;
        break;
      }
    }

    if(!done) {
      m_n.push_back(n);
      m_l.push_back(l);
      m_j.push_back(j);
      m_s.push_back(s);
    }  
  }

  else{
    m_n.push_back(n);
    m_l.push_back(l);
    m_j.push_back(j);
    m_s.push_back(s);
  }
}

////////////////////////////////////////////////////////////////////////////////
unsigned int ShellModelState::GetNumberOfOrbital(){
  return m_n.size();
}

////////////////////////////////////////////////////////////////////////////////
std::string ShellModelState::GetOrbitalString(unsigned int i,std::string option){
  std::ostringstream os;
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

  if(option=="")
    os << m_n[i] << shell << (m_j[i]*2) << "/2";

  else if(option=="latex")
    os << m_n[i] << shell <<"$_{" << (m_j[i]*2) << "/2}$";


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
  unsigned int orbital = 4294967295;
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
double ShellModelState::CumulativeSFDifference(ShellModelState& state){
  double diff=0;
  bool touch = false;
  unsigned int mysize1 = m_n.size();
  unsigned int mysize2 = state.GetNumberOfOrbital();

  for(unsigned int i = 0 ; i < mysize1 ; i++){
    for(unsigned int j = 0 ; j < mysize2 ; j++){
      if(m_n[i] == state.GetOrbitalN(j)
          && m_l[i] == state.GetOrbitalL(j)
          && m_j[i] == state.GetOrbitalJ(j)){
        diff+=fabs(m_s[i]-state.GetOrbitalS(j));
        touch = true;
      }
    }
  }

  if(!touch) 
    diff = 10000;

  return diff;

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
DifferentialCrossSection ShellModelState::GetTotalCS(){
  DifferentialCrossSection CS;
  unsigned int mysize = m_n.size();
  for(unsigned int i = 0 ; i < mysize ; i++){
    DifferentialCrossSection temp = GetOrbitalCS(i);
    CS +=temp;  
  } 

  return CS;
}
////////////////////////////////////////////////////////////////////////////////
std::vector<DifferentialCrossSection> ShellModelState::GetAllOrbitalCS(int limit){
  std::vector<DifferentialCrossSection> CS;
  unsigned int mysize = m_n.size();

  if(limit < 0){
    for(unsigned int i = 0 ; i < mysize ; i++){
      DifferentialCrossSection temp = GetOrbitalCS(i);
      CS.push_back(temp); 
    } 
  }

  else{
    // retain the limit main contributions
    // if the contribution have the same L but different J
    // Only the main J orbital is kept
    int mycount = 0 ;
    std::set<int> L;
    for(unsigned int i = 0 ; i < mysize && mycount < limit; i++){
      if(L.find(GetOrbitalL(i))==L.end()){
        L.insert(GetOrbitalL(i));
        DifferentialCrossSection temp = GetOrbitalCS(i);
        CS.push_back(temp); 
        mycount++;
      }
    } 
  }


  return CS;
}

////////////////////////////////////////////////////////////////////////////////
void ShellModelState::RemoveSecondaryLContribution(){
  std::set<int> Lfound;
  std::vector<unsigned int>::iterator it_n = m_n.begin();;
  std::vector<unsigned int>::iterator it_l = m_l.begin();;
  std::vector<double>::iterator it_j = m_j.begin();
  std::vector<double>::iterator it_s = m_s.begin();
  bool done = false;
  for(it_n = m_n.begin(); it_n!=m_n.end(); ){
    if(Lfound.find(*it_l)==Lfound.end()){
      Lfound.insert(*it_l);
      it_n++;it_l++;it_j++;it_s++;
    }
    else{
    it_n=  m_n.erase(it_n);
    it_l=  m_l.erase(it_l);
    it_j=  m_j.erase(it_j);
    it_s=  m_s.erase(it_s);
    }
  }

}


////////////////////////////////////////////////////////////////////////////////
DifferentialCrossSection ShellModelState::GetOrbitalCS(unsigned int i){
  double initial_spin = 2.5;
  double BeamEnergy =  5;
  double QValue = 3.35;
  double final_spin = m_J;

  // Test if the .rstemp folder exist
  // If yes do nothing
  // if no create it

  DIR *dir;
  if((dir = opendir(".rstemp"))==NULL){
    system("mkdir .rstemp");
    if((dir = opendir(".rstemp"))==NULL){
      std::cout << "RS++ fail to create temporary folder. Check your read/write permission for the current folder" << std::endl;
    }
  }
  closedir(dir);

  std::ofstream Front_Input(".rstemp/in.front");
  Front_Input << "jjj" << std::endl;
  Front_Input << "pipo" << std::endl;
  Front_Input << 2 << std::endl;
  Front_Input << BeamEnergy << std::endl;
  Front_Input << 25 << " " << 11 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << "0 0 0" << std::endl;
  Front_Input << m_l[i] << " " << m_j[i] << std::endl;
  Front_Input << m_n[i]-1 << std::endl;
  Front_Input << 2 << std::endl;

  // unbound case:
  if( QValue - m_E < 0 )
    Front_Input << 0.01 << std::endl;
  else
    Front_Input << QValue-m_E << std::endl; 

  Front_Input << 1 << std::endl;
  Front_Input << initial_spin << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << 5 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << final_spin << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << 2 << std::endl;
  Front_Input << 2 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << 1 << std::endl;
  Front_Input << 1.25 << " " << 0.65 << std::endl;
  Front_Input << 6.0 << std::endl;
  Front_Input << 0 << std::endl;
  Front_Input << 0 << std::endl;

  Front_Input.close() ;

  system("(cd .rstemp && exec FRONT < in.front &> /dev/null)"); 
  system("(cd .rstemp && exec echo tran.jjj | TWOFNR &> /dev/null)");

  DifferentialCrossSection CS;
  CS.LoadFromTWOFNR(".rstemp/24.jjj");
  CS.Scale(m_s[i]);
  return CS;
}

