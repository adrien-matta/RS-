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
      void LoadCollectionFromOxbash(unsigned int NumberOfState,std::string LPE, std::string LPF);
      void SetReferenceEnergy(double Energy);
      void SetGroundState(double J, int P , int order);
      void AddState(RS::ShellModelState state);
      unsigned int GetNumberOfState();
      RS::ShellModelState GetState(unsigned int i);
      void SetName(std::string name);
      void SelectStateByTotalSF(double threshold);
      void SelectStateByTotalCS(double threshold);
      void SelectStateByStrength(double threshold);
      void SelectStateByMainSF(double threshold);
      void SelectStateByMainCS(double threshold);
      void Print(int status = -1);
  };
}


#endif
