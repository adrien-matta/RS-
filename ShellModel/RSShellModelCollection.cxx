// STL
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
// RS++
#include "RSShellModelCollection.h"
using namespace RS;

////////////////////////////////////////////////////////////////////////////////
ShellModelCollection::ShellModelCollection(){
}

////////////////////////////////////////////////////////////////////////////////
ShellModelCollection::ShellModelCollection(std::string name){
  m_name = name;
}

////////////////////////////////////////////////////////////////////////////////
ShellModelCollection::~ShellModelCollection(){
}

////////////////////////////////////////////////////////////////////////////////
void ShellModelCollection::AddState(RS::ShellModelState state){
  m_collection.push_back(state);
  m_status.push_back(1);
}

////////////////////////////////////////////////////////////////////////////////
void ShellModelCollection::LoadCollectionFromOxbash(unsigned int NumberOfStates,std::string LPE, std::string LSF){
  /// Look how much state are already loaded:
  unsigned int offset = m_collection.size();


  //////////////
  // LPE PART //
  //////////////
  // Load the LPE File
  std::fstream levelfile(LPE);
  if(!levelfile.is_open()){
    std::cout << "Fail to open " << LPE << std::endl;
    exit(1);
  }

  double thespin ;
  double theparity;
  std::string LineBuffer,buffer;
  double gs = 0 ;
  // Read the file line by line
  while(getline(levelfile,LineBuffer)){
    // Indicate the part of the file that contain the spin of the states
    if(LineBuffer.find("j :")!=std::string::npos){
      std::istringstream os(LineBuffer);
      os >> buffer >> buffer >> thespin ;
    }

    // Indicate the part of the file that contain the spin of the states
    if(LineBuffer.find("+v")!=std::string::npos){
      theparity = 1 ;
    }

    if(LineBuffer.find("-v")!=std::string::npos){
      theparity = -1 ;
    }


    // Indicate the part of the file that contain the gs energy:
    if(LineBuffer.compare(0,16," gs energy     :")==0){
      std::istringstream os(LineBuffer);
      os >> buffer >> buffer >> buffer >> gs ;
    }

    // Indicate the part of the file that contain the states energies
    if(LineBuffer.compare(0,18," no   energy level")==0){
      // Take the next line
      getline(levelfile,LineBuffer);// empty
      getline(levelfile,LineBuffer);// empty
      int count = 0 ;
      while( getline(levelfile,LineBuffer) && count<NumberOfStates ){
        double E;
        std::istringstream os(LineBuffer);
        os >> buffer >> E ;
        ShellModelState state(E+gs,thespin,theparity,count+1);
        AddState(state);
        // Skip the next line
        if(LineBuffer.compare(0,18,"                   ")==0)
          getline(levelfile,LineBuffer);
        count++;
      }
    }
  }
  levelfile.close();


  //////////////
  // LSF PART //
  //////////////
  levelfile.open(LSF);
  if(!levelfile.is_open()){
    std::cout << "SF file " << LSF << " not found" << std::endl;
    return;
  }

  double theSF=0;
  double locSF=0;
  std::string linestart="";
  bool status = false;
  bool first = true;
  bool core  = false;
  int counter1 = 0 ;
  int counter2 = 0 ;
  std::vector<std::string> locOrbital;
  // Read the file line by line
  std::vector<int> vn;
  std::vector<int> vl;
  std::vector<int> vj;

  while(getline(levelfile,LineBuffer)){
    // Indicate the part of the file that describe the orbital
    if(LineBuffer.compare(0,18," -- core state ---")==0||(core&&LineBuffer.compare(0,5," 2j2t")!=0)){
      core = true ;
      std::string myLine = LineBuffer.substr(LineBuffer.find_first_of("("),LineBuffer.find_last_of(")")); 
      // search and replace all the parenthesis with space
      while(myLine.find("(")!=std::string::npos){
        myLine.replace(myLine.find("("),1," ");
      }

      while(myLine.find(")")!=std::string::npos){
        myLine.replace(myLine.find(")"),1," ");
      }


      std::istringstream los(myLine);
      double n,j,l,dump;


      // Create the list of orbital used in the file
      while(los >> n >> l >> j){
        vn.push_back(n);
        vl.push_back(l);
        vj.push_back(j);
      }
    }                           

    // Indicate the part of the file that contain the SF
    if(LineBuffer.compare(0,5," 2j2t")==0){
      core = false;
      // Take the next line
      getline(levelfile,LineBuffer);// gs energy, ignore
      //    getline(levelfile,LineBuffer);


      for(unsigned int i = 0 ; i < NumberOfStates ; i++){
        getline(levelfile,LineBuffer);
        if(LineBuffer.length()<vn.size()) getline(levelfile,LineBuffer); 

        if(LineBuffer.compare(0,6," total")==0){
          return;
        }

        std::istringstream os(LineBuffer);

        for(unsigned int j = 0 ; j < 10 ; j++)
          os >> buffer;

        for(unsigned j = 0 ; j < vn.size() ; j++){
          os >> theSF;
          if(theSF>0)
            m_collection[offset+i].AddOrbital(vn[j],vl[j],vj[j]/2.,theSF);
        }
      }
    }
  }
}

////////////////////////////////////////////////////////////////////////////////
unsigned int ShellModelCollection::GetNumberOfState(){
  return m_collection.size();
}

////////////////////////////////////////////////////////////////////////////////
void ShellModelCollection::SetName(std::string name){
  m_name = name;
}
////////////////////////////////////////////////////////////////////////////////
void ShellModelCollection::SelectStateByTotalSF(double threshold){
  unsigned int mysize = m_collection.size();
  for(unsigned int i = 0 ; i < mysize ; i++){
    double totalSF = 0;
    for(unsigned int l = 0 ; l < 8 ; l++){
      totalSF += m_collection[i].GetSumOfSForL(l);
    }

    if(totalSF > threshold){
      m_status[i]=1;
    }

    else
      m_status[i]=0;
  }
}
////////////////////////////////////////////////////////////////////////////////
void ShellModelCollection::SelectStateByTotalCS(double threshold){
  double initial_spin = 2.5;
  double BeamEnergy =  10;
  double QValue = 3.35; 
  unsigned int mysize = m_collection.size();
  for(unsigned int i = 0 ; i < mysize ; i++){
    double totalCS = 0;
    double final_spin = m_collection[i].GetJ();
    unsigned int nbrorb = m_collection[i].GetNumberOfOrbital();
    for(unsigned int orb = 0 ; orb < nbrorb ; orb++ ){
      std::ofstream Front_Input("in.front");
      Front_Input << "jjj" << std::endl;
      Front_Input << "pipo" << std::endl;
      Front_Input << 2 << std::endl;
      Front_Input << BeamEnergy << std::endl;
      Front_Input << 26 << " " << 11 << std::endl;
      Front_Input << 1 << std::endl;
      Front_Input << 1 << std::endl;
      Front_Input << "0 0 0" << std::endl;
      Front_Input << m_collection[i].GetOrbitalL(orb) << " " << m_collection[i].GetOrbitalJ(orb) << std::endl;
      Front_Input << m_collection[i].GetOrbitalN(orb)-1 << std::endl;
      Front_Input << 2 << std::endl;
      Front_Input << QValue - m_collection[i].GetEnergy() << std::endl;
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

      system("FRONT < in.front");
      system("echo tran.jjj | TWOFNR > dump.txt");

      double CS = 0 ;
      std::ifstream csfile("21.jjj");
      std::vector<double> x,y;
      double xx , yy, dump;
      
      while(csfile>> yy >> xx >> dump){
        CS+=sin(yy*M_PI/180.)*(2*M_PI)*(M_PI/180.)*xx;
      }

      totalCS +=CS*m_collection[i].GetOrbitalS(orb);
    }
    
    if(totalCS>threshold)
      m_status[i] = 1;
    else
      m_status[i] = 0;
  
  }

}

////////////////////////////////////////////////////////////////////////////////
void ShellModelCollection::SelectStateByStrength(double threshold){
  unsigned int mysize = m_collection.size();
  for(unsigned int i = 0 ; i < mysize ; i++){
    unsigned int NumberOfOrbital =  m_collection[i].GetNumberOfOrbital();
    double Strength = 0 ;
    for(unsigned int orb = 0 ; orb < NumberOfOrbital ; orb++){
      double factor = (2*m_collection[i].GetJ()+1)
        /(2*m_collection[i].GetOrbitalJ(orb)+1);
      Strength += m_collection[i].GetOrbitalS(orb)*factor;
    }

    if(Strength > threshold){
      m_status[i]=1;
    }

    else
      m_status[i]=0;
  }

}

////////////////////////////////////////////////////////////////////////////////
void ShellModelCollection::SelectStateByMainSF(double threshold){
  unsigned int mysize = m_collection.size();
  for(unsigned int i = 0 ; i < mysize ; i++){
    unsigned int Orbital =  m_collection[i].GetMainOrbital();

    if(Orbital == 1000)  m_status[i]=0;

    else if(m_collection[i].GetSumOfSForL( m_collection[i].GetOrbitalL(Orbital))> threshold){
      m_status[i]=1;
    }
    else
      m_status[i]=0;
  }
}
////////////////////////////////////////////////////////////////////////////////
void ShellModelCollection::SelectStateByMainCS(double threshold){
}
////////////////////////////////////////////////////////////////////////////////
void ShellModelCollection::Print(int status){
  // Print all
  if(status < 0 ){
    unsigned int mysize = m_collection.size();
    std::cout << "//////////////////// " << std::endl;
    std::cout << "Name : " << m_name << std::endl;
    std::cout << "Number of State : " << mysize << std::endl;

    for(unsigned int i = 0 ; i < mysize ; i++){
      m_collection[i].Print();
    }
    std::cout << "//////////////////// " << std::endl;
  }

  if(status >= 0 ){
    // Count number of selected state
    unsigned int Selected = 0;
    unsigned int mysize = m_collection.size();

    for(unsigned int i = 0 ; i < mysize ; i++)
      if(status == m_status[i])
        Selected++;


    std::cout << "//////////////////// " << std::endl;
    std::cout << "Name : " << m_name << std::endl;
    std::cout << "Number of State : " << mysize << std::endl;
    std::cout << "Number of Selected State: " << Selected << std::endl;

    for(unsigned int i = 0 ; i < mysize ; i++){
      if(status == m_status[i])
        m_collection[i].Print();
    }
  }
}
////////////////////////////////////////////////////////////////////////////////
void ShellModelCollection::SetReferenceEnergy(double Energy){
  if(Energy < 0 ) Energy = -Energy;

  unsigned int mysize = m_collection.size();
  for(unsigned int i = 0 ; i < mysize ; i++){
    m_collection[i].SetEnergy(m_collection[i].GetEnergy()+Energy);
  }

}
////////////////////////////////////////////////////////////////////////////////
void ShellModelCollection::SetGroundState(double J, int P , int order){
  unsigned int mysize = m_collection.size();
  for(unsigned int i = 0 ; i < mysize ; i++){
    if(m_collection[i].GetJ()==J
        && m_collection[i].GetParity()==P
        && m_collection[i].GetOrder()==order){
      SetReferenceEnergy( m_collection[i].GetEnergy());
    }
  }


}
////////////////////////////////////////////////////////////////////////////////