#include <chrono>
#include <thread>
#include <functional>
#include <memory>
#include <algorithm>
#include "boost/program_options.hpp"
#include "multiDomainWave.hpp"
#include "scalarWavePlots.hpp"


template <typename T>
struct checkHistoryEval
{
  template<typename U,void (U::*f)(const std::vector<double> &, double)> struct checker{};
  
  template<typename U>
  static std::true_type test(checker<U,&U::operator()>*);

  template<typename U>
  static std::false_type test(...);
  
  static const bool value = std::is_same<std::true_type, decltype(test<T>(nullptr))>::value;
};


template <typename T>
struct checkWaveEval
{
  template<typename U,void (U::*f)(const std::vector<double> &, std::vector<double> &, const double)> struct checker{};
  
  template<typename U>
  static std::true_type test(checker<U,&U::operator()>*);

  template<typename U>
  static std::false_type test(...);
  
  static const bool value = std::is_same<std::true_type, decltype(test<T>(nullptr))>::value;
};

  

void odeEvolve(std::vector<double> initial, auto &wave, double duration, double stepSize,
		     auto &waveHist){
  static_assert(checkHistoryEval<typename std::remove_reference<decltype(waveHist)>::type >::value,
		"odeEvolve was passed an invalid History with which to record");
  static_assert(checkWaveEval<typename std::remove_reference<decltype(wave)>::type >::value,
		"odeEvolve was passed an invalid wave to evolve");
  boost::numeric::odeint::runge_kutta4<std::vector<double>> rk;
  size_t steps =  boost::numeric::odeint::integrate_const(rk,*wave,initial,0.0,duration,stepSize,waveHist);
  printf("\ncompleted! number of steps: %d\n",(int)steps);
}

int main(int argv, char * args[])
{

  // Flags to add in:
  //  - initial data:
  //  -- gaussian
  //  -- sinusoid
  //  -- fast sinusoid
  //  - visualization
  //  -- animation (duration)
  //  -- mode plot (number of modes)
  //  - data output (spectral modes vs pseudospectral collocation points)
  //  -- print
  //  -- dump to file
  //  - number of domains
  //  - simulation params
  //  -- duration
  //  -- timestep
  //  -- spectral order

  // Initialize wave
  boost::program_options::options_description desc("Options:");
  desc.add_options()
    ("help","show this help message")
    ("id",boost::program_options::value<std::string>(),"specify initial data type (sin,fastsin,pulse)")
    ("bc",boost::program_options::value<std::string>(),"specify right boundary condition (transmit,reflect)")
    ("data","dump time-series collocation data to stdout")
    ("dom",boost::program_options::value<int>(),"specify number of domains")
    ("dur",boost::program_options::value<double>(),"duration of simulation")
    ("step",boost::program_options::value<double>(),"size of simulation timestep")
    ("type",boost::program_options::value<std::string>(),"type of spectral simulation (coll,dg)")
    ("ord",boost::program_options::value<int>(),"spectral order")
    ("no-vis","turn off default visualizations");

  boost::program_options::variables_map vars;
  boost::program_options::store(boost::program_options::parse_command_line(argv,args,desc),vars);

  if(vars.count("help")) {
    desc.print(std::cout);
    return 1;
  }

  // we can afford to only have a single variable function as we'll specify the ID to be right-going
  std::function<double(double)> boundData = [](double x){return cos(2*(x));};
  if(vars.count("id"))
    {
      if(vars["id"].as<std::string>() == "fastsin")
	boundData = [](double x){return cos(10*x);};
      else if(vars["id"].as<std::string>() == "pulse")
	boundData = [](double x){return pow(2,-5.0*pow(x,2));};
      else if(vars["id"].as<std::string>() != "sin")
	printf("id specified but does not match flags, defaulting to sin\n");
    }
  bool isReflecting = false;
  if(vars.count("bc"))
    {
      if(vars["bc"].as<std::string>() == "reflect")
	isReflecting = true;
      else if(vars["bc"].as<std::string>() != "transmit")
	printf("bc specified but does not match flags, defualting to transmit\n");	
    }
  bool isDG = true;
  if(vars.count("type"))
    {
      if(vars["type"].as<std::string>() == "coll")
	isReflecting = true;
      else if(vars["type"].as<std::string>() != "dg")
	printf("bc specified but does not match flags, defualting to transmit\n");
    }
  
  bool dumpData = (bool)(vars.count("data"));
  bool vis = !(bool)(vars.count("no-vis"));
  int doms;
  double duration;
  double step;
  int order;
  vars.count("dom") ? doms = vars["dom"].as<int>() : doms=2;
  vars.count("dur") ? duration = vars["dur"].as<double>() : duration = 10;
  vars.count("step") ? step = vars["step"].as<double>() : step = .01;
  vars.count("ord") ? order = vars["ord"].as<int>() : order = 20;
  
  std::vector<double> x = std::vector<double>();
  std::vector<double> times;
  times.push_back(0.0);

  std::vector<int> orders(doms,order);

  //TODO: Split all the options out to various functions
  std::vector<elementStateHistory> states;
	  
  std::vector<std::shared_ptr<std::vector<double>>> abscissas(doms);
  std::vector<std::shared_ptr<std::vector<double>>> weights(doms);
  std::vector<std::vector<double>> *derivmat;	  
  std::vector<std::shared_ptr<matrix<double>>> DMats(doms);
  std::transform(orders.begin(),orders.end(),abscissas.begin(),
		 isDG ? legendreTools::generateAbscissas : legendreTools::generateGLAbscissas);
  std::transform(orders.begin(),orders.end(),abscissas.begin(),weights.begin(),
		 isDG ? legendreTools::generateWeights : legendreTools::generateGLWeights);

  derivmat = legendreTools::generateDMat(order[d],abscissas[d],
					 legendreTools::generateBaryWeights(order[d],abscissas[d]));
  DMats.push_back(legendreTools::generateDMat(order[d],abscissas[d],
					      legendreTools::generateBaryWeights(order[d],abscissas[d])));

  /**/
  if(collocationMethod)
    {

      for(int d=0;d<numDomains;d++)
	{
	  abscissas.push_back(legendreTools::generateGLAbscissas(order[d]));
	  weights.push_back(legendreTools::generateGLWeights(order[d],abscissas[d]));
	}
      for(int d=0;d<numDomains;d++)
	{
	  for(int i=0;i<order[d];i++)
	    x.push_back(cos(-2.0*(abscissas[d]->at(i) - abscissas[d]->at(0) + 2.0*d)));
	  for(int i=0;i<order[d];i++)
	    x.push_back(-cos(-2.0*(abscissas[d]->at(i) - abscissas[d]->at(0) + 2.0*d)));
	}
      int elstart=0;
      for(int d=0;d<numDomains;d++)
	{
	  states.push_back(elementStateHistory());
	  for(int i=0;i<2;i++)
	    {
	      states[d].functionStates.push_back(functionStateHistory());
	      states[d].functionStates[i].timeStates.push_back(scalarFunction(order[i],abscissas[i],weights[i],DMats[i],
									      std::vector<double>(x.begin()+elstart
												  ,x.begin()+elstart+order[i])));
	      states[d].functionStates.push_back(functionStateHistory());
	      states[d].functionStates[i].timeStates.push_back(scalarFunction(order[i],abscissas[i],weights[i],DMats[i],
									      std::vector<double>(x.begin()+elstart+order[i]
												  ,x.begin()+elstart+2*order[i])));
	    }
	  elstart+=2*order[d];
	}
      multiStateHistory waveHist(order,numDomains,2,states,times);

	      
      collTransmittingMultiWave* wave = new collTransmittingMultiWave(order,abscissas,weights,DMats,stepSize,
								      numDomains,true);
      printf("initiating ode integrator\n");
	  
      odeEvolve(x,wave,10,.01,waveHist);

    }
  else
    {
	      
      for(int d=0;d<numDomains;d++)
	{
	  abscissas.push_back(legendreTools::generateAbscissas(order[d]));
	  weights.push_back(legendreTools::generateWeights(order[d],abscissas[d]));
	  dmatCatch = legendreTools::generateDMat(order[d],abscissas[d],
						  legendreTools::generateBaryWeights(order[d],abscissas[d]));
	  DMats.push_back(new matrix<double>(order[d],*dmatCatch));
	}
      for(int d=0;d<numDomains;d++)
	{
	  for(int i=0;i<order[d];i++)
	    x.push_back(cos(-2.0*(abscissas[d]->at(i) + 1.0  + 2.0*d)));
	  for(int i=0;i<order[d];i++)
	    x.push_back(-cos(-2.0*(abscissas[d]->at(i) + 1.0 + 2.0*d)));
	}
      int elstart=0;
      for(int d=0;d<numDomains;d++)
	{
	  states.push_back(elementStateHistory());
	  for(int i=0;i<2;i++)
	    {
	      states[d].functionStates.push_back(functionStateHistory());
	      states[d].functionStates[i].timeStates.push_back(
							       scalarFunction(order[d],abscissas[d],weights[d],DMats[d],
									      std::vector<double>(x.begin()+elstart+i*order[d]
												  ,x.begin()+elstart+(i+1)*order[d])));
	    }
	  elstart+=2*order[d];
	}
      multiStateHistory waveHist(order,numDomains,2,states,times);


	      
      DGTransmittingMultiWave* wave = new DGTransmittingMultiWave(order,abscissas,weights,DMats,stepSize,
								  numDomains,true);
      printf("initiating ode integrator\n");
      odeEvolve(x,wave,20,.01,waveHist);
    }
  
  //TODO: make this auto-tune to the parameters, make more general

  printf("Computing legendre modes (summing quadratures)...\n");
  std::vector<scalarFunction> plotAccumulator;
  for(int i = 0;i<states[0].functionStates[0].timeStates.size();i+=1)
    {
      plotAccumulator= std::vector<scalarFunction>();
      for(int j=0;j<states.size();j++)
	plotAccumulator.push_back(states[j].functionStates[0].timeStates[i]);
      multiPlotWaveandDeriv(plotAccumulator);
    }	  
}
  

