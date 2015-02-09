#ifndef __RSExperimentalCrossSection__
#define __RSExperimentalCrossSection__

// RS
#include "RSDifferentialCrossSection.h"
using namespace RS;

// STL
#include<string>
#include<vector>

// ROOT
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TH1.h"


namespace RS{

  class ExperimentalCrossSection{
   public:
      ExperimentalCrossSection();
      ~ExperimentalCrossSection();
 
   public: // Data Manipulation
      void LoadDataFromASCII(std::string);
      // Set Data from various root object
      void SetData(TGraph*);
      void SetData(TH1*);
      // Get Data
      TGraphAsymmErrors* GetData();

   private:
      TGraphAsymmErrors* m_CS;

   public: // Data Normalisation
      void Scale(double);
      void Divide(TH1*);
      void Divide(TF1*);
      void Divide(TGraph*);

   public: // Fit
      // Return the normalised CS as a new object
       DifferentialCrossSection FindNormalisation(std::vector<DifferentialCrossSection>&,std::vector<double>& Norm,std::vector<double>& errNorm);
       double ComputeChi2(DifferentialCrossSection&);
       double ComputeChi2Fit(const double* parameter);

   private:
       std::vector<DifferentialCrossSection> m_ThOriginal;
       std::vector<DifferentialCrossSection> m_ThScaled;

  };
}


#endif
