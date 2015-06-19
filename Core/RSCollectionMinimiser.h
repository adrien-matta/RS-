// STL
#include<vector>


// Root
#include"TGraph.h"
#include"TGraphErrors.h"

namespace RS{
  class CollectionMinimiser{
    public: 
      CollectionMinimiser();
      ~CollectionMinimiser();

    public:
      double Minimise(std::vector<TGraphErrors*> Collection, std::vector< std::vector<TGraph*> > Bases, std::vector<double>& param, std::vector<double>& err, double lmin=-1, double lmax=-1);
      double ComputeOverallChi2(const double* param);
      double ComputeOneChi2(unsigned int i, const double* param);

    private:
      std::vector<TGraphErrors*> m_Collection;
      std::vector< std::vector<TGraph*> > m_Bases;
      double m_min;
      double m_max;
  };
}
