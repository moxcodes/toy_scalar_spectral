#include <vector>
#include "math.h"
#include <boost/math/special_functions/legendre.hpp>
#include "legendreTools.hpp"
#include "matrix.hpp"

/* This class will define the behavior of a scalar function of the (x,t) variables
 *  - Represented as a spectral decomposition
 *  - Represented in a finite domain
 *
 *
 */

#ifndef SCALARFUNCTION
#define SCALARFUNCTION

class scalarFunction
{
public:
  std::vector<double> collocationData;
  // For convenience, each scalarFunction stores a pointer to the abscissas, weights, and first derivative matrix
  std::vector<double>* abscissas;
  std::vector<double>* weights;
  matrix<double>* DMat;
  // Lazily populated spectral data:
  std::vector<double>* spectralData;
  int n;

  scalarFunction(int order, std::vector<double>* inAbscissas, std::vector<double>* inWeights, matrix<double>* inDMat)
    : n(order), abscissas(inAbscissas), weights(inWeights), DMat(inDMat), spectralData(NULL) {}

  scalarFunction(int order, std::vector<double>* inAbscissas, std::vector<double>* inWeights, matrix<double>* inDMat,
		 std::vector<double> inCollocationData)
    : n(order), abscissas(inAbscissas), weights(inWeights), DMat(inDMat), spectralData(NULL), collocationData(inCollocationData) {}

  
  // for at collocation points
  double atCP(int i);

  double at(int i);

  double at(double x,std::vector<double>* baryWeights);
  
  double dxCP(int i);

  double dx(double x,std::vector<double>* baryWeights);
  
  double ddxCP(int i);

  // These each involve the quadrature sum, which is done 
  
  double at(double x);

  double dx(double x);

  double ddx(double x);

  void quadSum();

};
#endif
