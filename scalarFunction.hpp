#include <vector>
#include <memory>
#include "math.h"
#include <boost/math/special_functions/legendre.hpp>
#include "legendreTools.hpp"
#include "matrix.hpp"

#ifndef SCALARFUNCTION
#define SCALARFUNCTION

/// Single-variable spectral function

/// This class represents a single-variable function on a single Legendre
/// domain. It stores the collocation data, and has utilities to generate
/// spectral coefficients as well as compute the function value, first, and
/// second derivatives at collocation points and off.
class scalarFunction
{
public:
  std::vector<double> collocationData; ///< vector storing the list of collocation values, length n
  std::shared_ptr<std::vector<double>> abscissas; ///< a pointer to the abscissas of order n
  std::shared_ptr<std::vector<double>> weights; ///< a pointer to the weights of order n
  std::shared_ptr<matrix<double>> DMat; ///< a pointer to the derivative matrix for the abscissas used
  std::shared_ptr<std::vector<double>> spectralData; ///< The spectral coefficients for the function. Lazily populated, often null 
  int n; ///< The legendre order of the function


  /// Scalar function constructor. This takes as arguments the order of the
  /// function, the abscissas, the weights, and the derivative matrix
  /// pointers. This is to ensure that these objects are generated and duplicates
  /// are stored as little as possible
  /// \param order order of the Legendre expansion, should match length of
  /// abscissas, weights, derivative matrix
  /// \param inAbscissas a shared pointer to the vector of abscissas
  /// \param inWeights a shared pointer to the vector of quadrature weights (not
  /// interpolation weights)
  /// \param inDMat a shared pointer to a the derivative matrix associated with
  /// the also provided abscissas
  scalarFunction(int order, std::shared_ptr<std::vector<double>> inAbscissas, std::shared_ptr<std::vector<double>> inWeights,
		 std::shared_ptr<matrix<double>> inDMat)
    : n(order), abscissas(inAbscissas), weights(inWeights), DMat(inDMat), spectralData(NULL) {}

  /// Scalar function constructor. This takes as arguments the order of the
  /// function, the abscissas, the weights, and the derivative matrix
  /// pointers. This is to ensure that these objects are generated and
  /// duplicates are stored as little as possible. 
  /// \param order order of the Legendre expansion, should match length of
  /// abscissas, weights, derivative matrix
  /// \param inAbscissas a shared pointer to the vector of abscissas
  /// \param inWeights a shared pointer to the vector of quadrature weights (not
  /// interpolation weights)
  /// \param inDMat a shared pointer to a the derivative matrix associated with
  /// the also provided abscissas
  /// \param inCollocationData a vector of input collocation data to initialize
  scalarFunction(int order, std::shared_ptr<std::vector<double>> inAbscissas, std::shared_ptr<std::vector<double>> inWeights,
		 std::shared_ptr<matrix<double>> inDMat,std::vector<double> inCollocationData)
    : n(order), abscissas(inAbscissas), weights(inWeights), DMat(inDMat), spectralData(NULL), collocationData(inCollocationData) {}

  
  /// Evaluates the scalar function at a collocation point
  /// \param i the integer index of the collocation point to evaluate
  /// \return the value of the function at the point \f$f(x^n_i)\f$
  double at(int i);

  /// Evaluates the value of the scalar function at a collocation point. Uses
  /// closed-form Legendre function, so can have noise at +/-1.0 abscissas, and
  /// causes quadSum() to be run.
  /// \param x the value of the point to be evaluated
  /// \return the value of the function at the point \f$f(x^n_i)\f$
  double at(double x);
  
  /// Evaluates the scalar function at an arbitrary point using interpolation
  /// according to barycentric weights
  /// \param x the value of the point to be evaluated
  /// \param baryWeights a pointer to the vector of barycentric weights for interpolation
  /// \return the value of the function at the point \f$f(x)\f$
  double at(double x,std::shared_ptr<std::vector<double>> baryWeights);

  /// Evaluates the first derivative of the scalar function at a collocation
  /// point
  /// \param i the integer index of the collocation point to evaluate
  /// \return the value of the first derivative of the function at the point
  /// \f$f^\prime(x^n_i)\f$
  double dx(int i);

  /// Evaluates the first derivative of the scalar function at a collocation
  /// point. Uses closed-form Legendre function, so can have noise at +/-1.0
  /// abscissas, and causes quadSum() to be run.
  /// \param x the value of the point to be evaluated
  /// \return the value of the first derivative of the function at the point
  /// \f$f^\prime(x^n_i)\f$
  double dx(double x);

  /// Evaluates the first derivative of the scalar function at an arbitrary
  /// point using interpolation according to barycentric weights
  /// \param x the value of the point to be evaluated
  /// \param baryWeights a pointer to the vector of barycentric weights for interpolation
  /// \return the value of the first derivative of the function at the point \f$f^\prime(x)\f$  
  double dx(double x,std::shared_ptr<std::vector<double>> baryWeights);

  /// Evaluates the second derivative of the scalar function at a collocation
  /// point
  /// \param i the integer index of the collocation point to evaluate
  /// \return the value of the second derivative of the function at the point
  /// \f$f^{\prime \prime}(x^n_i)\f$  
  double ddx(int i);

  /// Evaluates the second derivative of the scalar function at a collocation
  /// point. Uses closed-form Legendre function, so can have noise at +/-1.0
  /// abscissas, and causes quadSum() to be run.
  /// \param x the value of the point to be evaluated
  /// \return the value of the second derivative of the function at the point
  /// \f$f^{\prime \prime}(x^n_i)\f$  
  double ddx(double x);

  /// Performs the quadrature sum to populate the spectralData using the
  /// abscissas, weights, and collocationData. Called by functions which
  /// evaluate the function at arbitrary points.
  void quadSum();

};
#endif
