#ifndef __RSShellModelState__
#define __RSShellModelState__
#include<string>
#include<vector>

// RS
#include"RSDifferentialCrossSection.h"

namespace RS{

  class ShellModelState{
    public:
      ShellModelState();
      ShellModelState(double E, double J, int P,int order);
      ~ShellModelState();

    private:
      double m_E;
      double m_J;
      int m_P;
      int m_order;
      std::vector<unsigned int> m_n;
      std::vector<unsigned int> m_l;
      std::vector<double> m_j;
      std::vector<double> m_s;

    public:
      void AddOrbital(unsigned int n, unsigned int l, double j, double s);
      double GetEnergy();
      void SetEnergy(double Energy);
      double GetJ();
      int GetParity();
      int GetOrder();
      unsigned int GetNumberOfOrbital();
      std::string GetOrbitalString(unsigned int i,std::string option="");
      unsigned int GetOrbitalN(unsigned int i);
      unsigned int GetOrbitalL(unsigned int i);
      double GetOrbitalJ(unsigned int i);
      double GetOrbitalS(unsigned int i);
     
      DifferentialCrossSection GetTotalCS();
      DifferentialCrossSection GetOrbitalCS(unsigned int i);
      std::vector<DifferentialCrossSection> GetAllOrbitalCS(int limit =-1);
      
      unsigned int GetMainOrbital();
      double GetSumOfSForL(unsigned int l);
      // Compute the difference in SF for two given state
      // Used to Match state between collection
      double CumulativeSFDifference(ShellModelState&);
      
      // Loop over the orbital and for each L keep the main orbital
      void RemoveSecondaryLContribution();

      void Print();
  };
}
#endif
