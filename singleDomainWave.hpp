#include <vector>
#include <boost/math/tools/roots.hpp>
#include <boost/numeric/odeint.hpp>
#include "scalarFunction.hpp"
#include <stdio.h>


//TODO: turn this into a scalarfunction list type
struct stateHistory
{
  std::vector<std::vector<scalarFunction>>& m_states;
  std::vector<double >& m_times;
  int n;
  
  stateHistory(int length, std::vector<std::vector<scalarFunction>> &states , std::vector<double > &times )
    : n(length), m_states( states ) , m_times( times ) { }

  void operator()( const std::vector<double> &x ,double t )
  {
    int numScalars = x.size()/n;
    for(int i=0;i<numScalars;i++)
      {
	scalarFunction newVal = m_states[i].back();
	newVal.collocationData = std::vector<double>(x.begin() + i*n, x.begin() + i*n + n);
	m_states[i].push_back(newVal);
      }
    m_times.push_back(t);
  }
};


class singleDomainWave{
public:
  int n;
  std::vector<double>* abscissas;
  std::vector<double>* weights;
  matrix<double>* DMat;
  bool verbose;
  
  singleDomainWave(int ord,   std::vector<double>* in_abscissas,std::vector<double>* in_weights,
		   matrix<double> * in_DMat, double timestep, bool in_verbose)
  : n(ord),abscissas(in_abscissas),weights(in_weights),DMat(in_DMat), verbose(in_verbose) {}
  
  virtual void operator () (const std::vector<double> &x, std::vector<double> &dxdt, const double t) {}
};




// In this, we have a 2d state vector
// x[0-(n-1)] -  pi
// x[n-(2n-1)] - psi

class transmittingWave: public singleDomainWave{
  
public:
  transmittingWave(int ord,   std::vector<double>* in_abscissas,std::vector<double>* in_weights,
		   matrix<double> * in_DMat, double timestep, bool in_verbose)
    : singleDomainWave(ord,in_abscissas,in_weights,in_DMat,timestep,in_verbose) {}


  void operator() ( const std::vector<double> &x , std::vector<double> &dxdt , const double t)
  {
       
    std::vector<double>* pi = new std::vector<double>(x.begin(),x.begin()+n);
    std::vector<double>* psi = new std::vector<double>(x.begin()+n,x.begin()+ 2* n);

    dxdt[0] = -2*sin(2*(t));
    dxdt[n] = 2*sin(2*(t));
    
    std::vector<double> dpsi = ((*DMat)*(*psi));
    std::vector<double> dpi = ((*DMat)*(*pi));
    for(int i=1;i<n-1;i++)
      {
	dxdt[i] = dpsi[i];
	dxdt[n+i] = dpi[i];
      }

    dxdt[n-1]=-dpi[n-1];
    dxdt[2*n - 1] =-dpsi[n-1];

    if((int)t==t &&  verbose)
      printf("t=%f...",t);
    
  }
};

class advection: public singleDomainWave{
  
public:
  advection(int ord,   std::vector<double>* in_abscissas,std::vector<double>* in_weights,
	    matrix<double> * in_DMat, double timestep, bool in_verbose)
    : singleDomainWave(ord,in_abscissas,in_weights,in_DMat,timestep,in_verbose) {}

  void operator() ( const std::vector<double> &x , std::vector<double> &dxdt , const double t)
  {    
    std::vector<double>* pi = new std::vector<double>(x.begin(),x.begin()+n);

    dxdt[0] = -2*sin(2*(t));
    
    std::vector<double> dpi = ((*DMat)*(*pi));
    for(int i=1;i<n;i++)
	dxdt[i] = -dpi[i];
   
    if((int)t == t && verbose)
      printf("t=%f...",t);
  }
};

