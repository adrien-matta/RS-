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
      double Integrate(double low=0, double up=180);
      void Scale(const double& scale);
      DifferentialCrossSection& operator+=(DifferentialCrossSection& left);
      DifferentialCrossSection& operator-=(DifferentialCrossSection& left);

    public: // ROOT
      double Eval(const double& x) ;
      TGraph* GetTGraph();
      bool m_modified;  
    private: // ROOT
    TGraph* m_g ;

  };
  inline DifferentialCrossSection operator+(DifferentialCrossSection& left, DifferentialCrossSection& right);
  inline DifferentialCrossSection operator-(DifferentialCrossSection& left, DifferentialCrossSection& right);


}

#endif
