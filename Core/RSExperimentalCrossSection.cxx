#include"RSExperimentalCrossSection.h"
using namespace RS;

// STL
#include<iostream>
#include<fstream>
#include<sstream>
// ROOT
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
////////////////////////////////////////////////////////////////////////////////
ExperimentalCrossSection::ExperimentalCrossSection(){
  m_CS=0;
}
////////////////////////////////////////////////////////////////////////////////
 ExperimentalCrossSection::~ExperimentalCrossSection(){
}
////////////////////////////////////////////////////////////////////////////////
void ExperimentalCrossSection::LoadDataFromASCII(std::string filename){
  std::ifstream file(filename.c_str()); 
  if(!file.is_open()){
    std::cout << "ERROR : File " << filename << " not found" << std::endl;
    exit(1);
  }

  // First analyse the file to check what kind of error are presents:
  std::string LineBuffer;
  getline(file,LineBuffer);
  std::istringstream line(LineBuffer);
  double buffer;
  unsigned int mycount = 0 ;
  while(line >> buffer)
     mycount++;

  std::vector<double> x,exh,exl,y,eyh,eyl;
  double bx,bexh,bexl,by,beyh,beyl;

  file.close();
  file.open(filename.c_str());

  if(mycount == 2){
    while(file >> bx >> by){
      x.push_back(bx);
      exh.push_back(0);
      exl.push_back(0);
      y.push_back(by);
      eyh.push_back(0);
      eyl.push_back(0);
    }
  }
  
  else if(mycount == 3){
    while(file >> bx >> by >> beyl){
      x.push_back(bx);
      exh.push_back(0);
      exl.push_back(0);
      y.push_back(by);
      eyh.push_back(beyl);
      eyl.push_back(beyl);
    }
  }

  else if(mycount == 4){
    while(file >> bx >> by >> beyl >> beyh){
      x.push_back(bx);
      exh.push_back(0);
      exl.push_back(0);
      y.push_back(by);
      eyh.push_back(beyl);
      eyl.push_back(beyh);
    }
  }

  else if(mycount == 5){
    while(file >> bx >> bexl >> by >> beyl >> beyh){
      x.push_back(bx);
      exh.push_back(bexl);
      exl.push_back(bexl);
      y.push_back(by);
      eyh.push_back(beyl);
      eyl.push_back(beyh);
    }
  }

  else if(mycount == 6){
    while(file >> bx >> bexl >> bexh >> by >> beyl >> beyh){
      x.push_back(bx);
      exh.push_back(bexl);
      exl.push_back(bexh);
      y.push_back(by);
      eyh.push_back(beyl);
      eyl.push_back(beyh);
    }
  }

  file.close();
  if(m_CS)
    delete m_CS;
  m_CS = new TGraphAsymmErrors(x.size(), &x[0], &y[0],&exl[0],&exh[0], &eyl[0],&eyh[0]);
}
////////////////////////////////////////////////////////////////////////////////
void ExperimentalCrossSection::SetData(TGraph*){
}
////////////////////////////////////////////////////////////////////////////////
void ExperimentalCrossSection::SetData(TH1*){
}
////////////////////////////////////////////////////////////////////////////////
TGraphAsymmErrors* ExperimentalCrossSection::GetData(){
  return m_CS;
}
////////////////////////////////////////////////////////////////////////////////
void ExperimentalCrossSection::Scale(double){
}
////////////////////////////////////////////////////////////////////////////////
void ExperimentalCrossSection::Divide(TH1*){
}
////////////////////////////////////////////////////////////////////////////////
void ExperimentalCrossSection::Divide(TF1*){
}
////////////////////////////////////////////////////////////////////////////////
void ExperimentalCrossSection::Divide(TGraph*){
}
////////////////////////////////////////////////////////////////////////////////
DifferentialCrossSection ExperimentalCrossSection::FindNormalisation(std::vector<DifferentialCrossSection> Th,std::vector<double>& Norm,std::vector<double>& errNorm){
  // Make a local copy of the function
  m_ThOriginal.clear();
  m_ThOriginal = Th;
  unsigned int mysize = Th.size();
  Norm.resize(mysize,0.5);
  errNorm.resize(mysize,0.5);

  const char* minName ="Minuit";const char* algoName="Fumili2";
  ROOT::Math::Minimizer* min =
    ROOT::Math::Factory::CreateMinimizer(minName, algoName);
  min->SetValidError(true);

  // create funciton wrapper for minimizer
  // a IMultiGenFunction type
   ROOT::Math::Functor f(this,&ExperimentalCrossSection::ComputeChi2Fit, mysize);
  min->SetFunction(f);
  
  double* parameter = new double[mysize];
  for(unsigned int i = 0 ; i < mysize ; i++){
    parameter[i] = Norm[i];
    char name[4];
    sprintf(name,"S%d",i+1);
    min->SetLimitedVariable(i,name,0.5,0.1,0,1);
  }
  
  // do the minimization
  min->Minimize();
  const double *xs = min->X();
  const double *err = min->Errors(); 
 // Return the Fitted CS
  DifferentialCrossSection cs;
  for(unsigned int i = 0 ; i < mysize ; i++){
    Norm[i] = xs[i];
    errNorm[i] = err[i]; 
    Th[i].Scale(xs[i]);
    cs+=Th[i];
  }

  // Clean up
  m_ThOriginal.clear();
  m_ThScaled.clear(); 

  return cs;

}
////////////////////////////////////////////////////////////////////////////////
double ExperimentalCrossSection::ComputeChi2(DifferentialCrossSection& DCS){
  TGraph* th = DCS.GetTGraph();
  double chi2 = 0;
  double chi ;
  int mysize = m_CS->GetN();
  double* x = m_CS->GetX();
  double* y = m_CS->GetY();
  for(int i = 0 ; i < mysize ; i++){
    chi =  (y[i]-th->Eval(x[i]))/m_CS->GetErrorY(i);
    chi2+=chi*chi;
  }

  return chi2; 
}
////////////////////////////////////////////////////////////////////////////////
double ExperimentalCrossSection::ComputeChi2Fit(const double* parameter){
  unsigned int mysize = m_ThOriginal.size();
  m_ThScaled.resize(mysize);
  DifferentialCrossSection CS;
    for(unsigned int i = 0 ; i < mysize ; i++){
     m_ThScaled[i]= m_ThOriginal[i];
     m_ThScaled[i].Scale(parameter[i]);
     CS= CS + m_ThScaled[i];
  }
  return ComputeChi2(CS); 
}
