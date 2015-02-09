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
      void LoadFromTWOFNR(std::string filename);
      double Integrate(double low=-1, double up=200);
      void Scale(const double& scale);
      DifferentialCrossSection& operator+=(DifferentialCrossSection&);
      DifferentialCrossSection& operator-=(DifferentialCrossSection& );

    public: // ROOT
      double Eval(const double& x);
      TGraph* GetTGraph();
      bool m_modified; 

    private: // ROOT
      TGraph* m_g ;

  };
  DifferentialCrossSection operator+(DifferentialCrossSection& , DifferentialCrossSection& );
  DifferentialCrossSection operator-(DifferentialCrossSection& , DifferentialCrossSection& );


}

#endif
