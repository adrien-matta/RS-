#ifndef __RSShellModelCollection__
#define __RSShellModelCollection__

// STL
#include<string>
#include<vector>

// RS++
#include "RSShellModelState.h"

namespace RS{
  class ShellModelCollection{
    public:
      ShellModelCollection();
      ShellModelCollection(std::string name);
      ~ShellModelCollection();

    private:
      std::vector<RS::ShellModelState> m_collection;
      std::vector<unsigned int> m_status;
      std::string m_name;

    public:
      std::vector< std::vector <unsigned int> > MatchCollection(ShellModelCollection& Collection,double limit, unsigned int status);
      void LoadCollectionFromOxbash(unsigned int NumberOfState,std::string LPE, std::string LPF);
      void LoadCollectionFromSimpleFile(std::string FileName);
      void SetReferenceEnergy(double Energy);
      void SetGroundState(double J, int P , int order);
      void AddState(RS::ShellModelState state);
      unsigned int GetStatus(unsigned int i);
      unsigned int GetNumberOfState();
      RS::ShellModelState GetState(unsigned int i);
      void SetName(std::string name);
      void SelectStateByTotalSF(double threshold);
      void SelectStateByTotalCS(double threshold);
      void SelectStateByStrength(double threshold);
      void SelectStateByMainSF(double threshold);
      void SelectStateByMainCS(double threshold);
      void SelectStateByParity(int Parity);
      void Print(int status = -1);
  };
}


#endif
