#ifndef __RSDifferentialCrossSection__
#define __RSDifferentialCrossSection__
// STL
#include<string>
#include<vector>


// ROOT
#include"TGraph.h"

namespace RS{
  class DifferentialCrossSection{
    public:
      DifferentialCrossSection();
      ~DifferentialCrossSection();

    private:
      std::vector<double> m_x;
      std::vector<double> m_y;

    public:
      void LoadFromASCII(std::string filename);
      double Integrate(double low=-1, double up=-1);
      void Scale(double scale);
      DifferentialCrossSection& operator+=(const DifferentialCrossSection& left);
    public: // ROOT
      TGraph* GetTGraph();
  
    private: // ROOT
    TGraph* g ;

  };
  inline DifferentialCrossSection operator+(DifferentialCrossSection& left, const DifferentialCrossSection& right);
  inline DifferentialCrossSection operator-(DifferentialCrossSection& left, const DifferentialCrossSection& right);


}

#endif
