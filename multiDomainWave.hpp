#include <vector>
#include <boost/math/tools/roots.hpp>
#include <boost/numeric/odeint.hpp>
#include "scalarFunction.hpp"
#include <stdio.h>

/*
There's a bit of a row-column ordering here. Innermost structure is the time series, so longest.
ascii diagram:
(func)
||||   ||||   ||||        ||||
----   ----   ----        ----
 |      |      |           |
 --------      -------------
 (element)           |
     |               |
     -----------------
     (history)
*/
// TODO: rethink this structuring. At least overload [] in each
struct functionStateHistory
{
  std::vector<scalarFunction> timeStates;
  functionStateHistory(){}
};

struct elementStateHistory
{
  std::vector<functionStateHistory> functionStates;
  elementStateHistory(){}
};


struct multiStateHistory
{
  std::vector<elementStateHistory> &elementStates;
  std::vector<double>& m_times;
  std::vector<int> n;
  int doms;
  int funcs;

  //assumes full tree is intantiated with an initial data element
  multiStateHistory(std::vector<int> in_n, int domains, int functions, std::vector<elementStateHistory> &states, std::vector<double> &times)
    : n(in_n), elementStates(states), doms(domains), funcs(functions), m_times(times){}

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
  int test(double i){return 0;}
};

class multiDomainWave{
public:
  std::vector<int> n;
  int doms;
  std::vector<std::vector<double>*> abscissas;
  std::vector<std::vector<double>*> weights;
  std::vector<matrix<double>*> DMats;
  bool verbose;

  multiDomainWave(std::vector<int> ord, std::vector<std::vector<double>*> in_abscissas,
		  std::vector<std::vector<double>*> in_weights, std::vector<matrix<double>*> in_DMats,
		  double timestep, int domains, bool in_verbose)
    : n(ord), abscissas(in_abscissas),weights(in_weights),DMats(in_DMats), verbose(in_verbose), doms(domains) {}

  virtual void operator () (const std::vector<double> &x, std::vector<double> &dxdt, const double t){}
};

class collTransmittingMultiWave: public multiDomainWave{

public:
  collTransmittingMultiWave(std::vector<int> ord, std::vector<std::vector<double>*> in_abscissas,
			    std::vector<std::vector<double>*> in_weights,std::vector<matrix<double>*> in_DMats,
			    double timestep, int domains, bool in_verbose)
    : multiDomainWave(ord,in_abscissas,in_weights,in_DMats,timestep,domains,in_verbose) {}
  
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
    leftderivspi[0] = -rightderivspi[0] -2*2.0*sin(2.0*t);
    leftderivspsi[0] = -rightderivspsi[0] + 2*2.0*sin(2.0*t);
    rightderivspi[doms] = - leftderivspi[doms] - 2*leftderivspsi[doms];
    rightderivspsi[doms] = -leftderivspsi[doms] - 2*leftderivspi[doms];
    //edges evolve via averaging
    elstart=0;
    for(int i=0;i<doms;i++)
      {
	dxdt[elstart] = (leftderivspi[i] + rightderivspi[i])/2.0;
	dxdt[elstart + n[i] - 1] = (leftderivspi[i+1] + rightderivspi[i+1])/2.0;
	dxdt[elstart + n[i]] = (leftderivspsi[i] + rightderivspsi[i])/2.0;
	dxdt[elstart + 2*n[i] - 1] = (leftderivspsi[i+1] + rightderivspsi[i+1])/2.0;
	elstart+=2*n[i];
      }
  }

  std::vector<double> bulkEvolve(const std::vector<double> &x, std::vector<double> &dxdt,int el,int elstart)
  {
    std::vector<double>* pi = new std::vector<double>(x.begin()+elstart,x.begin()+elstart+n[el]);
    std::vector<double>* psi = new std::vector<double>(x.begin()+elstart+n[el],x.begin()+elstart+2*n[el]);
    std::vector<double> dpsi = ((*DMats[el])*(*psi));
    std::vector<double> dpi = ((*DMats[el])*(*pi));
    for(int i=1;i<n[el]-1;i++)
      {
	dxdt[elstart + i] = dpsi[i];
	dxdt[elstart + n[el]+i] = dpi[i];
      }
    return std::vector<double>({dpsi[0],dpi[0],dpsi[n[el]-1],dpi[n[el]-1]});
  }
};


class DGTransmittingMultiWave: public multiDomainWave{
public:
  // required for evaluating the flux values
  std::vector<std::vector<double>> leftInterpolant;
  std::vector<std::vector<double>> rightInterpolant;
  std::vector<std::vector<double>*> baryWeights;

  std::vector<matrix<double>*> DMatsHat;
  
  DGTransmittingMultiWave(std::vector<int> ord, std::vector<std::vector<double>*> in_abscissas,
		  std::vector<std::vector<double>*> in_weights,std::vector<matrix<double>*> in_DMats,
			double timestep, int domains, bool in_verbose)
    : multiDomainWave(ord,in_abscissas,in_weights,in_DMats,timestep,domains,in_verbose)
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

    rightfluxpi.push_back(cos(2.0*t));
    rightfluxpsi.push_back(-cos(2.0*t));

    for(int d=0;d<doms;d++)
      {
	leftfluxpi.push_back((pis[d].at(-1.0) + psis[d].at(-1.0))/2.0);
     	rightfluxpi.push_back((pis[d].at(1.0) - psis[d].at(1.0))/2.0);	
	leftfluxpsi.push_back((pis[d].at(-1.0) + psis[d].at(-1.0))/2.0);
	rightfluxpsi.push_back((-pis[d].at(1.0) + psis[d].at(1.0))/2.0);
      }

    
    leftfluxpi.push_back(0);
    leftfluxpsi.push_back(0);

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
  }

};
