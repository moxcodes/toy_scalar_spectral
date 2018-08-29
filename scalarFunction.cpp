#include "scalarFunction.hpp"
#include <stdio.h>


double scalarFunction::at(int i){
  return collocationData[i];}

double scalarFunction::dx(int i){
  return ((*DMat)*collocationData)[i];}

double scalarFunction::ddx(int i){
  return ((*DMat)*(*DMat)*collocationData)[i];}


double scalarFunction::at(double x){
  if(spectralData == NULL)
    quadSum();
  double val=0;
  for(int i=0;i<n;i++)
      val+= spectralData->at(i)* boost::math::legendre_p(i,x);
  return val;}


double scalarFunction::dx(double x){
  if(spectralData == NULL)
    quadSum();
  double val=0;
  for(int i=0;i<n;i++)
    val+= spectralData->at(i)*legendreTools::legendreDeriv(i,x);
  return val;
}

double scalarFunction::at(double x,std::shared_ptr<std::vector<double>> baryWeights){
  double den = 0;
  double num = 0;
  double prod;
  for(int i=0;i<n;i++)
    {
      prod = baryWeights->at(i)/(x - abscissas->at(i));
      num += prod*(collocationData[i]);
      den += prod;
    }
  return num/den;
}




double scalarFunction::dx(double x,std::shared_ptr<std::vector<double>> baryWeights){
  double den = 0;
  double num = 0;
  double pointval = at(x);
  double prod;
  for(int i=0;i<n;i++)
    {
      prod = baryWeights->at(i)/(x - abscissas->at(i));
      num += prod*(pointval - collocationData[i])/(x - abscissas->at(i));
      den += prod;
    }
  return num/den;
}



double scalarFunction::ddx(double x)
{
  if(spectralData == NULL)
    quadSum();
  double val=0;
  for(int i=0;i<n;i++)
    val+= spectralData->at(i)*legendreTools::legendreDDeriv(i,x);
  return val;
}


void scalarFunction::quadSum()
{
  if(spectralData == NULL)
    spectralData = std::shared_ptr<std::vector<double>>(new std::vector<double>((size_t)n));
  
  for(int i=0;i<n;i++)
    {
      double legi = 0;
      for(int j=0;j<n;j++)
	  legi+=weights->at(j)*collocationData[j]*boost::math::legendre_p(i,abscissas->at(j));
      spectralData->at(i) = legi*(2*i + 1)/2;
    }
  return;
}
