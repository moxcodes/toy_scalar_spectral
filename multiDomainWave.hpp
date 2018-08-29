#include <vector>
#include <boost/math/tools/roots.hpp>
#include <boost/numeric/odeint.hpp>
#include "scalarFunction.hpp"
#include <stdio.h>

/// A structure for holding a function state history, which just hold a vector
/// of scalarFunction representing the time series for a single function on a
/// single domain.
struct functionStateHistory
{
  std::vector<scalarFunction> timeStates; ///< the time series data for a single function 
  functionStateHistory(){}
};

/// A structure for holding the history of a single domain, which just holds a
/// vector of function state histories (as each domain can have multiple
/// functions)
struct elementStateHistory
{
  std::vector<functionStateHistory> functionStates;///< a vector of function histories
  elementStateHistory(){}
};


/// A structure for storing the history of a full wave evolution, which can have
/// multiple functions in each domain, and potentially several domains. This
/// structure can be used as the history object to be passed to boost ode
/// libraries.
struct multiStateHistory
{
  std::vector<elementStateHistory> &elementStates; ///< the vector of histories for the several domains
  std::vector<double>& m_times; ///< a vector of times, which should align with the states at the innermost of the domain heirarchies
  std::vector<int> n; ///< legendre order, so number of collocation points in each function
  int doms; ///< number of domains
  int funcs; ///< number of functions

  /// history constructor. Importantly, the elementStateHistory MUST be populated
  /// with the initial data for the routine to work properly, as scalarFunctions
  /// store more data than can be provided by the boost ode library
  /// \param in_n Legendre order of approximation
  /// \param domains number of domains in the simulation
  /// \param functions number of functions in each domain  
  /// \param states the Histories, already populated with initial data
  /// scalarFunctions - these will be populated by reference
  /// \param times the times, populated with initial time - these will be
  /// populated by reference
  multiStateHistory(std::vector<int> in_n, int domains, int functions,
		    std::vector<elementStateHistory> &states, std::vector<double> &times)
    : n(in_n), elementStates(states), doms(domains), funcs(functions), m_times(times){}


  /// Storage operator for use with boost ode integrators. Takes a flat input
  /// and organizes it into the spectral data hierarchy, stores it in the
  /// history values. Organization is assumed to follow structure of (domain 0,
  /// function 0);(domain 0, function 1);(domain 1, function 0);(domain 1,
  /// function 1)...
  /// \param x raw flattened ode data
  /// \param t simulation time
  void operator()(const std::vector<double> &x, double t)
  {
    int elstart =0;
    for(int d=0;d<doms;d++)
      {
	for(int f=0;f<funcs;f++)
	  {
	    scalarFunction newVal = elementStates[d].functionStates[f].timeStates.back();
	    newVal.collocationData = std::vector<double>(x.begin()+elstart+f*n[d],x.begin()+elstart+(f+1)*n[d]);
	    elementStates[d].functionStates[f].timeStates.push_back(newVal);
	  }
	elstart+=2*n[d];
      }
    m_times.push_back(t);
  }
};


/// Parent class for the various wave function implementations
class multiDomainWave{
public:
  std::vector<int> n;///< spectral order of the evolution
  int doms; ///< number of domains
  std::vector<std::shared_ptr<std::vector<double>>> abscissas; ///< a vector of abscissas storage, one for each domain
  std::vector<std::shared_ptr<std::vector<double>>> weights; ///< a vector of weight storage, one for each domain
  std::vector<std::shared_ptr<matrix<double>>> DMats; ///< a vector of derivative matrix storage, one for each domain
  std::function<double(double)> boundData; ///< a function for the left boundary data
  bool verbose;///< a flag for outputting status checkpoints to stdout

  /// The generic wave constructor, takes in much data about wave options
  /// \param ord Legendre order of simulation
  /// \param in_abscissas vector of abscissa data, one for each domain
  /// \param in_weights vector of weight data, one for each domain
  /// \param in_DMats vector of derivative matrix data, one for each domain
  /// \param domains number of domains
  /// \param in_boundData a function representing information about the driving boundary
  /// \param in_verbose whether the function should output status checkpounts
  multiDomainWave(std::vector<int> ord, std::vector<std::shared_ptr<std::vector<double>>> in_abscissas,
		  std::vector<std::shared_ptr<std::vector<double>>> in_weights,
		  std::vector<std::shared_ptr<matrix<double>>> in_DMats,
		  int domains, std::function<double(double)> in_boundData, bool in_verbose)
    : n(ord), abscissas(in_abscissas),weights(in_weights),DMats(in_DMats), verbose(in_verbose), boundData(in_boundData), doms(domains) {}

  /// virtual evolution operator - to be overwritten in all inherited classes
  virtual void operator () (const std::vector<double> &x, std::vector<double> &dxdt, const double t){}
};


/// Continuous boundary multiple-domain wave simulation class

/// A specialization class from multiDomainWave for the collocation point wave
/// multiple-domain simulation. This multi-domain method is intended for use
/// with Gauss-Lobatto abscissas, and evolves by ensuring consistency between
/// the -1 and 1 abscissas at neighboring domains.
class collTransmittingMultiWave: public multiDomainWave{

public:
  bool reflect; ///< true if reflecting right bound, false if transmitting

  /// The collocation wave constructor, takes in much data about wave options
  /// \param ord Legendre order of simulation
  /// \param in_abscissas vector of abscissa data, one for each domain
  /// \param in_weights vector of weight data, one for each domain
  /// \param in_DMats vector of derivative matrix data, one for each domain
  /// \param domains number of domains
  /// \param in_boundData a function representing the first time derivative used for left bound
  /// \param isReflecting true if right bound should reflect, false if transmit
  /// \param in_verbose whether the function should output status checkpounts
  collTransmittingMultiWave(std::vector<int> ord, std::vector<std::shared_ptr<std::vector<double>>> in_abscissas,
			    std::vector<std::shared_ptr<std::vector<double>>> in_weights,
			    std::vector<std::shared_ptr<matrix<double>>> in_DMats, int domains,
			    std::function<double(double)> in_boundData, bool isReflecting, bool in_verbose)
    : multiDomainWave(ord,in_abscissas,in_weights,in_DMats,domains,in_boundData,in_verbose), reflect(isReflecting) {}

  /// Wave evolution operator, for use in boost ode libraries. This gives the
  /// first derivative of each collocation point with respect to time by
  /// computing the first derivatives in each domain, and using neighboring
  /// domains to infer the derivative at the shared boundary points, which are
  /// constrained to evolve identically.
  /// \param x the set of flattened collocation points
  /// \param dxdt the set of first derivatives with respect to time - populated
  /// by this function as return parameter
  /// \param t simulation time of the timestep considered
  /// \sa bulkEvolve
  void operator() ( const std::vector<double> &x, std::vector<double> &dxdt, const double t)
  {
    // 2 steps : first evolve the bulk of each domain and get out the presumed time dependence
    //           at each interface, then impose agreement at each interface by averaging.
    //           This is slightly different from the Continuous Galerkin method, but only
    //           in the details of where the weight factors appear.

    
    // bulk evolve and get boundary vals
    std::vector<double> rightderivspi(doms+1);
    std::vector<double> rightderivspsi(doms+1);
    std::vector<double> leftderivspi(doms+1);
    std::vector<double> leftderivspsi(doms+1);
    //elstart tracks the position in the array of the start of the element, as they are of nonuniform size
    std::vector<double> derivsCatcher;
    int elstart = 0;
    for(int i=0;i<doms;i++)
      {
	derivsCatcher = bulkEvolve(x,dxdt,i,elstart);
	rightderivspi[i]=derivsCatcher[0];
	rightderivspsi[i]=derivsCatcher[1];
	leftderivspi[i+1]=derivsCatcher[2];
	leftderivspsi[i+1]=derivsCatcher[3];
	elstart+=2*n[i];
      }
    //leftmost and rightmost parts need to be modified to obey boundary conditions
    leftderivspi[0] = -rightderivspi[0] + 2*boundData(t+1.0);
    leftderivspsi[0] = rightderivspi[0];
    rightderivspi[doms] = - leftderivspi[doms] - (reflect ? 0 : 2*leftderivspsi[doms]);
    rightderivspsi[doms] = leftderivspsi[doms];
    //edges evolve via averaging
    elstart=0;
    for(int i=0;i<doms;i++)
      {
	dxdt[elstart] = (leftderivspi[i] + rightderivspi[i])/2.0;
	dxdt[elstart + n[i]] = (leftderivspsi[i] + rightderivspsi[i])/2.0;
	dxdt[elstart + n[i] - 1] = (leftderivspi[i+1] + rightderivspi[i+1])/2.0;
	dxdt[elstart + 2*n[i] - 1] = (leftderivspsi[i+1] + rightderivspsi[i+1])/2.0;
	elstart+=2*n[i];
      }
    if((int)(t) == t && verbose)
	printf("simulation time t=%f\n",t);
  }


  /// Function for evolving the bulk of the individual elements
  /// \param x the full flattened collocation points (for all domains)
  /// \param dxdt the full set of collocation first derivative to be populated,
  /// partially populated (for a single domain) as a return parameter of this
  /// function
  /// \param el number of the domain element to be populated
  /// \param elstart the index in the full flattened data of the start of the
  /// element to be populated
  /// \return a vector of the boundary derivatives (left psi, left pi, right
  /// psi, right pi)
  std::vector<double> bulkEvolve(const std::vector<double> &x, std::vector<double> &dxdt,int el,int elstart)
  {
    std::shared_ptr<std::vector<double>> pi(new std::vector<double>(x.begin()+elstart,x.begin()+elstart+n[el]));
    std::shared_ptr<std::vector<double>> psi(new std::vector<double>(x.begin()+elstart+n[el],x.begin()+elstart+2*n[el]));
    std::vector<double> dpsi = ((*DMats[el])*(*psi));
    std::vector<double> dpi = ((*DMats[el])*(*pi));
    for(int i=0;i<n[el];i++)
      {
	dxdt[elstart + i] = dpsi[i];
	dxdt[elstart + n[el] + i] = dpi[i];
      }
    return std::vector<double>({dpsi[0],dpi[0],dpsi[n[el]-1],dpi[n[el]-1]});
  }
};


/// DG simulation class

/// A specialization class from multiDomainWave for the Discontinuous Galerkin
/// wave multiple-domain simulation. This multi-domain method is intended for
/// use with Legendre Gauss abscissas, and evolves by imposing numerical fluxes
/// between domains
class DGTransmittingMultiWave: public multiDomainWave{
public:
  std::vector<std::vector<double>> leftInterpolant;///< the interpolant values for the point -1 for each domain
  std::vector<std::vector<double>> rightInterpolant;///< the interpolant values for the point 1 for each domain
  std::vector<std::shared_ptr<std::vector<double>>> baryWeights;///< the barycentric weights data for each domain

  std::vector<matrix<double>*> DMatsHat; ///< the adjusted matrix datas for the DG computation
  bool reflect;///< true if right boundary should reflect, false if transmit


  /// The discontinous Galerkin wave constructor, takes in much data about wave
  /// options. In addition to the typical data, it requires the interpolation
  /// functions evaluated at +/-1 and the adjusted derivative matrix according
  /// to ratios of weights. 
  /// \param ord Legendre order of simulation
  /// \param in_abscissas vector of abscissa data, one for each domain
  /// \param in_weights vector of weight data, one for each domain
  /// \param in_DMats vector of derivative matrix data, one for each domain
  /// \param domains number of domains
  /// \param in_boundData a function representing the first time derivative used for left bound
  /// \param isReflecting true if right bound should reflect, false if transmit
  /// \param in_verbose whether the function should output status checkpounts  
  DGTransmittingMultiWave(std::vector<int> ord, std::vector<std::shared_ptr<std::vector<double>>> in_abscissas,
			  std::vector<std::shared_ptr<std::vector<double>>> in_weights,
			  std::vector<std::shared_ptr<matrix<double>>> in_DMats, int domains,
			  std::function<double(double)> in_boundData, bool isReflecting, bool in_verbose)
    : multiDomainWave(ord,in_abscissas,in_weights,in_DMats,domains,in_boundData,in_verbose),reflect(isReflecting)
  {
    double lefts;
    double leftt;
    double rights;
    double rightt;
    for(int d=0;d<doms;d++)
      {
	leftInterpolant.push_back(std::vector<double>());
	rightInterpolant.push_back(std::vector<double>());
	lefts=0;
	rights=0;
	baryWeights.push_back(legendreTools::generateBaryWeights(n[d],abscissas[d]));
	for(int i=0;i<n[d];i++)
	  {
	    leftt=baryWeights[d]->at(i)/(-1.0 - abscissas[d]->at(i));
	    rightt=baryWeights[d]->at(i)/(1.0 - abscissas[d]->at(i));
	    leftInterpolant[d].push_back(leftt);
	    rightInterpolant[d].push_back(rightt);
	    lefts+=leftt;
	    rights+=rightt;
	  }
	for(int i=0;i<n[d];i++)
	  {
	    leftInterpolant[d][i]=leftInterpolant[d][i]/lefts;
	    rightInterpolant[d][i]=rightInterpolant[d][i]/rights;
	  }
	DMatsHat.push_back(new matrix<double>(n[d],DMats[d]->matData));
	for(int i=0;i<n[d];i++)
	  for(int j=0;j<n[d];j++)
	    DMatsHat[d]->matData[i][j] = -DMats[d]->matData[j][i] *weights[d]->at(j)/weights[d]->at(i);
      }
  }

  /// Wave evolution operator, for use in boost ode libraries. This gives the
  /// first derivative of each collocation point with respect to time by
  /// computing the left and right fluxes at each domain boundary, then using
  /// those in the DG formulas found in [Kopriva ch 8]
  /// \param x the set of flattened collocation points
  /// \param dxdt the set of first derivatives with respect to time - populated
  /// by this function as return parameter
  /// \param t simulation time of the timestep considered
  void operator() ( const std::vector<double> &x, std::vector<double> &dxdt, const double t)
  {

    // 2 steps : first derive the left and right fluxes for each domain and variable, then use them to
    //           evolve each of the collocation points.

    std::vector<double> rightfluxpi;
    std::vector<double> rightfluxpsi;
    std::vector<double> leftfluxpi;
    std::vector<double> leftfluxpsi;

    std::vector<scalarFunction> pis;
    std::vector<scalarFunction> psis;

    int elstart=0;
    //generate the vectors of field values
    for(int d=0;d<doms;d++)
      {
	pis.push_back(scalarFunction(n[d],abscissas[d],weights[d],DMats[d],
				     std::vector<double>(x.begin()+ elstart,x.begin()+elstart+n[d])));
	psis.push_back(scalarFunction(n[d],abscissas[d],weights[d],DMats[d],
				      std::vector<double>(x.begin()+elstart+n[d],x.begin()+elstart+ 2*n[d])));
	elstart+=2*n[d];
      }

    rightfluxpi.push_back(boundData(t+1.0));
    rightfluxpsi.push_back(-boundData(t+1.0));

    for(int d=0;d<doms;d++)
      {
	leftfluxpi.push_back((pis[d].at(-1.0) + psis[d].at(-1.0))/2.0);
     	rightfluxpi.push_back((pis[d].at(1.0) - psis[d].at(1.0))/2.0);	
	leftfluxpsi.push_back((pis[d].at(-1.0) + psis[d].at(-1.0))/2.0);
	rightfluxpsi.push_back((-pis[d].at(1.0) + psis[d].at(1.0))/2.0);
      }

    
    leftfluxpi.push_back(reflect ? -rightfluxpi[doms]: 0);
    leftfluxpsi.push_back(reflect ? -rightfluxpi[doms]: 0);

    //step 2: evolve the bulk using flux vals
    std::vector<double> dpsi;
    std::vector<double> dpi;
    elstart=0;
    for(int d=0;d<doms;d++)
      {
	dpsi = ((*DMatsHat[d])*(psis[d].collocationData));
	dpi = ((*DMatsHat[d])*(pis[d].collocationData));
	for(int i=0;i<n[d];i++)
	  {

	    dxdt[elstart + i] = (dpsi[i]
				 + (leftfluxpi[d+1] - rightfluxpi[d+1])*rightInterpolant[d][i]/weights[d]->at(i)
				 -(leftfluxpi[d] - rightfluxpi[d])*leftInterpolant[d][i]/weights[d]->at(i));
	    dxdt[elstart + i + n[d]] = (dpi[i]
					+ (leftfluxpsi[d+1] - rightfluxpsi[d+1])*rightInterpolant[d][i]/weights[d]->at(i)
					-(leftfluxpsi[d] - rightfluxpsi[d])*leftInterpolant[d][i]/weights[d]->at(i)); 
	  }
	elstart+=2*n[d];
      }
  
    if((int)(t) == t && verbose)
      printf("simulation time t=%f\n",t);
  }
};
