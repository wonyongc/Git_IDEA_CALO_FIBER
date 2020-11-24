#include "VectorSmallestInterval.hh"
#include "TMath.h"

double FindSmallestInterval(double mean, double meanErr, double min, double max, std::vector<double>* vals, const double fraction, const bool verbosity)
{
  if( verbosity )
    std::cout << ">>>>>> FindSmallestInterval" << std::endl;


  std::sort(vals->begin(),vals->end());

  unsigned int nPoints = vals->size();
  unsigned int maxPoints = (unsigned int)(fraction * nPoints);

//   unsigned int minPoint = 0;
//   unsigned int maxPoint = 0;
  double delta = 999999.;
  for(unsigned int point = 0; point < nPoints-maxPoints; ++point)
  {
    double tmpMin = vals->at(point);
    double tmpMax = vals->at(point+maxPoints-1);
    if( tmpMax-tmpMin < delta )
    {
      delta = tmpMax - tmpMin;
      min = tmpMin;
      max = tmpMax;
    }
  }

  
  return delta;
  
}

double FindSmallestInterval(std::vector<double>* vals, const double fraction, const bool verbosity)
{
  if( verbosity )
    std::cout << ">>>>>> FindSmallestInterval" << std::endl;


  std::sort(vals->begin(),vals->end());

  unsigned int nPoints = vals->size();
  unsigned int maxPoints = (unsigned int)(fraction * nPoints);

  double delta = 9999999.;
  for(unsigned int point = 0; point < nPoints-maxPoints; ++point)
  {
    double tmpMin = vals->at(point);
    double tmpMax = vals->at(point+maxPoints-1);
    if( tmpMax-tmpMin < delta )
    {
      delta = tmpMax - tmpMin;
    }
  }

  
  return delta;
  
}



double poissonf(double*x, double*par)                                         
{                                                    
//   Double_t res=0.;
  Double_t xx=x[0];
  if (xx<=0) return  0;

  // Poisson distribution
  // par[1] - distribution parameter
  return par[0]*TMath::Power(par[1],xx)/TMath::Gamma(xx+1)/TMath::Exp(par[1]);

//   return par[0]*TMath::Poisson(x[0],par[1]);
}