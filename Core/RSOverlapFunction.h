#ifndef _OVERLAP_
#define _OVERLAP_
// STL
#include<vector>
#include<string>

// Root
#include"TGraph.h"

namespace RS{
  class OverlapFunction{
    public:
      OverlapFunction();
      OverlapFunction(std::string filename);
      ~OverlapFunction();

    private:
      std::vector<double> m_OverlapFunction_R;
      std::vector<double> m_OverlapFunction_OL;
      std::vector<double> m_OverlapFunction_E;
      TGraph* m_Graph;
      bool m_modified;
    public:
      inline std::vector<double> GetR(){return m_OverlapFunction_R;};
      inline std::vector<double> GetOL() {return m_OverlapFunction_OL;};
      inline std::vector<double> GetE() {return m_OverlapFunction_E;};

    public:
      void LoadFromFile(std::string filename);
      void Sort();
      void Normalise(double norm);
      void Multiply(TF1*);
      void Divide(TGraph*);
      void Add(TGraph*,double f=1);
      void Scale(double factor);
      double Integrate(double min=-1, double max=-1);
      TGraph* GetTGraph();
      double Eval(double R);
      OverlapFunction& operator+=(OverlapFunction&);
      OverlapFunction& operator-=(OverlapFunction& );

  
    public: // Fit
      RS::OverlapFunction FitWithOverlap(std::vector<RS::OverlapFunction>,std::vector<double>& param,std::vector<double>& err,double lmin=-1,double lmax=-1,std::vector<RS::OverlapFunction> Constant=std::vector<RS::OverlapFunction>());
      double ComputeChi2(OverlapFunction&);
      double ComputeChi2Fit(const double* param);

    private:
      std::vector<RS::OverlapFunction> m_ThOriginal;
      std::vector<RS::OverlapFunction> m_ThCurrent;
      std::vector<RS::OverlapFunction> m_ThConstant;

      double m_min;
      double m_max;
  
  };
}
#endif
