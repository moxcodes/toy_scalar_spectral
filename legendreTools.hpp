#include <vector>
#include <memory>
#include <algorithm>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/tools/roots.hpp>
#include <math.h>
#include "matrix.hpp"

#ifndef LEGENDRETOOLS_H
#define LEGENDRETOOLS_H

const double PI = 3.1415926535897932384626433832795028841971693993751058209749445923078164062;

const int PREC = 50;

//! Tools for Legendre polynomials, and corresponding spectral quantities
namespace legendreTools{

  /// Value of the first derivative of a particular Legendre polynomial order at a particular point
  /// \param n the order of Legendre polynomial to evaluate
  /// \param x the point at which to evaluate the derivative
  /// \return the first derivative, evaluated from the analytic formula associated with legendre polys
  ///         \f$P_n^\prime(x) = n \frac{x P_n(x) - P_{n-1}(x)}{x^2 - 1}\f$
  /// \sa legendreDDeriv
  static double legendreDeriv(const int n, const double x)
  {
    if(n==0)
      return 0;
    else
      return n*(x* boost::math::legendre_p(n,x) - boost::math::legendre_p(n-1,x))/(pow(x,2) - 1);
  }

  /// Value of the second derivative of a particular Legendre polynomial order at a particular point
  /// \param n the order of Legendre polynomial to evaluate
  /// \param x the point at which to evaluate the second derivative
  /// \return the first derivative, evaluated from the analytic formula associated with legendre polys
  ///         \f$P_n^{\prime\prime}(x) = n\frac{((n - 1) x^2 - n - 1)P_n(x) + 2 x P_{n - 1}(x)}{(x^2 - 1)^2}\f$
  /// \sa legendreDeriv
  static double legendreDDeriv(int n, double x)
  {
    return n*((((n - 1)* pow(x,2) - n - 1) *  boost::math::legendre_p(n,x) + 2 * x * boost::math::legendre_p(n-1,x))
	       /pow(pow(x,2) - 1,2));
  }


  
  /// function object for Legendre polynomial
  
  /// A function object representing the Legendre polynomial of a particular order
  /// with its derivative. Used primarily for Newton-Raphson.
  class legendre{
  public:
    int n; ///< order of the Legendre polynomial \f$P_n\f$

    /// Legendre polynomial constructor; stores the order of the polynomial
    legendre(const int order)
    {n=order;};

    /// Legendre polynomial and derivative evaluator
    /// \param x the point at which to evaluate the polynomial
    /// \return A two-element tuple of doubles, first element is the value of the polynomial,
    ///         second element is the value of the first derivative \f$(P_n(x),P^\prime_n(x))\f$
    std::tuple<double,double> operator()(const double x){
      return std::make_tuple<double,double>(boost::math::legendre_p(n,x),
					    legendreDeriv(n,x));}

    /// Value of the Legendre polynomial at a point
    /// \param x the point at which to evaluate the polynomial
    /// \return the value of \f$P_n(x)\f$
    double at(const double x){
      return boost::math::legendre_p(n,x);}

    /// Value of the first derivative of the Legendre polynomial at a point
    /// \param x the point at which to evaluate the derivative
    /// \return the value of \f$P^\prime_n(x)\f$
    /// \sa legendreDeriv
    double dx(const double x){
      return legendreDeriv(n,x);
    }
     
  };

  /// Function object for difference of legendre polynomial n and n-2
  
  /// A function object representing the function \f$q_n(x) \equiv P_n(x) - P_{n-2}(x) \f$
  /// with its derivative. Used primarily for Newton-Raphson in Gauss-Lobatto generation.
  class q{
  public:
    int n; ///< order of the function \f$q_n\f$

    /// Constructor, simply stores the order of the polynomial
    q(int order)
    {n=order;};

    /// \f$q_n\f$ polynomial and derivative evaluator
    /// \param x the point at which to evaluate the polynomial
    /// \return A two-element tuple of doubles, first element is the value of the polynomial,
    ///         second element is the value of the first derivative \f$(q_n(x),q^\prime_n(x))\f$
    std::tuple<double,double> operator()(double x){
      return std::make_tuple<double,double>(boost::math::legendre_p(n,x) - boost::math::legendre_p(n-2,x),
						      legendreDeriv(n,x) -legendreDeriv(n-2,x));}

    /// Value of the \f$q_n\f$ polynomial at a point
    /// \param x the point at which to evaluate the polynomial
    /// \return the value of \f$q_n(x)\f$
    double at(double x){
      return boost::math::legendre_p(n,x) - boost::math::legendre_p(n-2,x);}


    /// Value of the first derivative of the \f$q_n\f$ polynomial at a point
    /// \param x the point at which to evaluate the derivative
    /// \return the value of \f$q^\prime_n(x)\f$
    /// \sa legendreDeriv    
    double dx(double x){
      return legendreDeriv(n,x) -legendreDeriv(n-2,x);}
  };



  
  /// Generates the Gauss-Legendre Abscissas at a particular order in the
  /// approximation The abscissas determine the collocation points for a
  /// spectral approximation to obtain desired spectral convergence in
  /// evaluating pseudospectral evolution.  The abscissas are located at the
  /// zeros of [Kopriva; Numerical Recipes] \f$P_n(x)\f$, where n is the number of
  /// spectral points
  /// \param order the order of Legendre polynomials for which to generate abscissas
  /// \return a shared_ptr to the vector of abscissas \f$x_n\f$
  static std::shared_ptr<std::vector<double>> generateAbscissas(int order){
    std::shared_ptr<std::vector<double>> abscissas(new std::vector<double>());

    // generate a very rough approximation for the first; working backward
    // for Newton-Raphson nextAb < ab < prevAb
    double prevAb = 1;
    double ab= ((1 - 1/(double)pow(order,2) + 1/(double)pow(order,3))
		     * cos(PI*((double)(4 - 1))/((double)(4*order + 2))));
    double nextAb;
    for(int i = 1;i<=(order/2);i++)
      {	
	nextAb= ((1 - 1/pow(order,2) + 1/pow(order,3))
			  * cos(PI*(4*(i+1) - 1)/(4*order + 2)));
	legendre legendreN = legendre(order);
	abscissas->insert(abscissas->begin(),
			  boost::math::tools::newton_raphson_iterate(legendreN,ab,
								     nextAb,prevAb,PREC));
	prevAb=ab;
	ab=nextAb;
      }

    // half generated, other half directly inferred
    if(order%2 !=0)
      abscissas->insert(abscissas->begin(),0);
    for(int i = 1;i<=(order/2);i++)
	abscissas->insert(abscissas->begin()+i-1,-1*abscissas->at(abscissas->size()-i));
    return abscissas;
  };


  /// Generates the Gauss-Lobatto Abscissas at a particular order in the
  /// approximation The abscissas determine the collocation points for a
  /// spectral approximation to obtain desired spectral convergence in
  /// evaluating pseudospectral evolution.  The abscissas are located at the
  /// zeros of [Kopriva] \f$P_n(x) - P_{n-2}(x)\f$, where n is the number of spectral points
  /// \param order the order of Legendre polynomials for which to generate abscissas
  /// \return a shared_ptr to the vector of abscissas \f$x_i^n\f$
  static std::shared_ptr<std::vector<double>> generateGLAbscissas(int order){
    std::shared_ptr<std::vector<double>> abscissas(new std::vector<double>());
    // generate a very rough approximation for the first couple.
    // working forward, prevAb < ab < nextAb
    double prevAb= -1;
    double ab = - cos( (double)(5.0/4.0)*PI/((double) (order - 1.0))
		       - 3.0/((double)(8.0*(order - 1.0)*PI*( 5.0/4.0))));
    double nextAb;
    abscissas->push_back(-1);
    for(int j=1;j<(order)/2;j++)
      {
	nextAb= - cos( (double)(j+1 + 1.0/4.0)*PI/((double) (order - 1.0 ))
		       - 3.0/((double)(8.0*(order - 1.0)*PI*(j+1 + 1.0/4.0))));
        q qN = q(order);
	abscissas->push_back(boost::math::tools::newton_raphson_iterate(qN,ab,
									prevAb,nextAb,PREC));
	prevAb=ab;
	ab=nextAb;
      }

    // half generated, other half directly inferred
    if(order%2 !=0)
      abscissas->push_back(0.0);
    for(int i = 1;i<=(order/2);i++)
      abscissas->push_back(-abscissas->at((order/2) - i));
    return abscissas;
  };

  /// Generates the Gauss-Legendre Weights at a particular order in the
  /// approximation. The weights are coefficients in the quadrature sum
  /// to convert from the data at collocation points to Legendre polynomial coefficients.
  /// Each weight is given by [Kopriva] : \f$w^n_i = \frac{2}{(1 - (x^n_i)^2)*P_n(x_i^n)^2}\f$
  /// \param order the order of Legendre polynomials for which to generate weights for the same order
  /// \param abscissas a shared_ptr to the vector of Gauss-Legendre abscissas 
  /// \return a shared_ptr to the vector of weights \f$w_i^n\f$
  static std::shared_ptr<std::vector<double>> generateWeights(int order, std::shared_ptr<std::vector<double>> abscissas)
  {
    std::shared_ptr<std::vector<double>> weights(new std::vector<double>());
    for(int i=0;i<order;i++)
      weights->push_back(2.0/((1.0 - pow(abscissas->at(i),2))*pow(legendreDeriv(order,abscissas->at(i)),2)));
    return weights;
  }

  /// Generates the Gauss-Lobatto Weights at a particular order in the
  /// approximation. The weights are coefficients in the quadrature sum
  /// to convert from the data at collocation points to Legendre polynomial coefficients.
  /// Each weight is given by [Kopriva] : \f$w^n_i = \frac{2}{n(n-1)(P_{n-1}(x^n_i))^2}\f$
  /// \param order the order of Legendre polynomials for which to generate weights
  /// \param abscissas a shared_ptr to the vector of Gauss-Lobatto abscissas for the same order 
  /// \return a shared_ptr to the vector of weights \f$w_i^n\f$  
  static std::shared_ptr<std::vector<double>> generateGLWeights(int order, std::shared_ptr<std::vector<double>> abscissas)
  {
    std::shared_ptr<std::vector<double>> weights(new std::vector<double>());
    weights->push_back(2.0/((double)(order*(order - 1))));
    for(int i=1;i<order-1;i++)
      weights->push_back(2.0/((order*(order -1.0))*pow(boost::math::legendre_p(order-1,abscissas->at(i)),2)));
    weights->push_back(2.0/((double)(order*(order - 1))));
    return weights;
  }

  /// Generates the Barycentric Weights at a particular order in the
  /// approximation. The weights are coefficients in the sum to convert between
  /// collocation data and interpolant values in the interval.
  /// Each weight is given by [Kopriva] : \f$w^{B\,n}_i = \frac{1}{\prod_{i=0;i\ne j}^n (x_i^n - x_j^n)}\f$
  /// \param order the order of Legendre polynomials for which to generate weights
  /// \param abscissas a shared_ptr to the vector of abscissas of desired type for the same order
  /// \return a shared_ptr to the vector of abscissas \f$w_i^{B\,n}\f$  
  static std::shared_ptr<std::vector<double>> generateBaryWeights(int order, std::shared_ptr<std::vector<double>> abscissas)
  {
    std::shared_ptr<std::vector<double>> weights( new std::vector<double>(order,1.0));
    for(int i=0;i<order;i++)
      {
	for(int j=0;j<order;j++)
	  {
	    if(i!=j)
	      weights->at(i) = weights->at(i)*(abscissas->at(i) - abscissas->at(j));
	  }
      }
    for(int i=0;i<order;i++)
      {
	weights->at(i) = 1.0/weights->at(i);
      }
    return weights;
  }

  /// Generates the matrix with which the derivatives at collocation points can
  /// be computed (either for Gauss-Legendre or Gauss-Lobatto collocation). The
  /// matrix values are computed as per [Kopriva]: \f$D_{i j}^n =
  /// \frac{w^n_j}{w^n_i} \frac{1}{x^n_i - x^n_j}\f$ for \f$i \ne j\f$ and
  /// \f$D_{i i} = -\sum_{j=0;j\ne i}^n D_{i j}\f$
  /// \param order the order at which the approximation is evaluated
  /// \param abscissas a vector of abscissas for which the derivative matrix will be used
  /// \param bWeights the barycentric weights for the set of abscissas
  /// \return a shared_ptr to the matrix of abscissas
  static std::shared_ptr<matrix<double>> generateDMat(int order, std::shared_ptr<std::vector<double>> abscissas,
						      std::shared_ptr<std::vector<double>> bWeights)
  {
    std::vector<std::vector<double>> DMat;
    for(int i=0;i<order;i++)
	DMat.push_back(std::vector<double>(order));
    
    for(int i=0;i<order;i++)
      {
	DMat[i][i] = 0;
	for(int j=0;j<order;j++)
	  {
	    if(i!=j)
	      {
		DMat[i][j] = ((bWeights->at(j)/bWeights->at(i)) *(1.0/(abscissas->at(i) - abscissas->at(j))));
		DMat[i][i] = (DMat[i][i] - DMat[i][j]);
	      }
	  }
      }
    return std::shared_ptr<matrix<double>>(new matrix<double>(order,DMat));
  }
}

#endif
