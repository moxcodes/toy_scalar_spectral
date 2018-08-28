#include "singleDomainWave.hpp"


//TODO: fix specialization to boundary conditions

singleDomainWave::singleDomainWave(int order, double timestep, long double dur){
  n = order;
  duration=dur;
  abscissas = legendreTools::generateGLAbscissas(n);
  weights = legendreTools::generateGLWeights(n,abscissas);
  DMat = legendreTools::generateDMat(n,abscissas,legendreTools::generateBaryWeights(n,abscissas));
  currentTimestep=0;
  timestepSize=timestep;


  std::vector<double>* piQuad = new std::vector<long double>();
  for(int i=0;i<n;i++)
    piQuad->push_back(cos(-abscissas->at(i)));

  times = std::vector<double>();
  scalarFunction pi = scalarFunction(order);
  pi.quadSum(weights,abscissas,piQuad);
  piHist.push_back(pi);
  sleep(1);
}






void singleDomainWave::odeEvolve()
{
  std::vector<double> x = std::vector<long double>();
  for(int i=0;i<n;i++)
    x.push_back(cos(-2*(abscissas->at(i) - abscissas->at(0))));
  //used if we use the transmittingwave:
  for(int i=0;i<n;i++)
    x.push_back(-cos(-2*(abscissas->at(i) - abscissas->at(0))));
    
  std::vector<std::vector<double>> states=std::vector<std::vector<long double>>();
  transmittingWave tw(n,abscissas,weights,DMat,timestepSize);
  boost::numeric::odeint::runge_kutta4< std::vector<double> > rk;
  size_t steps = boost::numeric::odeint::integrate_const(rk,tw,x,0.0,duration,timestepSize,
							 stateHistory(states,singleDomainWave::times));
  printf("number of steps: %d\n",steps);
  printf("size of array: %d\n",states.size());
  std::vector<double>::const_iterator first;
  std::vector<double>::const_iterator last;
  std::vector<double>* vec;
  scalarFunction func(n);
  //load data back into phis and pis:

  for(int i=1;i<=steps;i++)
    {
      first=states[i].begin();
      last=states[i].begin()+n;
      vec=new std::vector<double>(first,last);
      func.quadSum(weights,abscissas,vec);     
      piHist.push_back(func);
    }
  return;
}



