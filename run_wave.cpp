#include <chrono>
#include <thread>
#include <functional>
#include <memory>
#include <algorithm>
#include <boost/program_options.hpp>
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
  size_t steps =  boost::numeric::odeint::integrate_const(rk,wave,initial,0.0,duration,stepSize,waveHist);
  printf("\ncompleted! number of steps: %d\n",(int)steps);
}

int main(int argv, char * args[])
{

  //initialize command line options
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
    ("no-vis","turn off default visualizations")
    ("verbose","turn on periodic status updates during simulation");

  boost::program_options::variables_map vars;
  boost::program_options::store(boost::program_options::parse_command_line(argv,args,desc),vars);

  if(vars.count("help")) {
    desc.print(std::cout);
    return 1;
  }

  //assign command line options to local variables for use 
  
  // we can afford to only have a single variable function as we'll specify the ID to be right-going
  std::function<double(double)> boundData = [](double x){return cos(2*(x));};
  std::function<double(double)> boundDatadx = [](double x){return - 2.0*sin(2.0*x);};
  if(vars.count("id"))
    {
      if(vars["id"].as<std::string>() == "fastsin")
	{
	  boundData = [](double x){return cos(10*x);};
	  boundDatadx = [](double x){return -10.0*cos(10*x);};
	}
      else if(vars["id"].as<std::string>() == "pulse")
	{
	  boundData = [](double x){return pow(2,-5.0*pow(x,2));};
	  boundDatadx = [](double x){return -10.0*x*pow(2,-5.0*pow(x,2));};
	}
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
        isDG = false;
      else if(vars["type"].as<std::string>() != "dg")
	printf("type specified but does not match flags, defualting to transmit\n");
    }
  
  bool dumpData = (bool)(vars.count("data"));
  bool verb = (bool)(vars.count("verbose"));
  bool vis = !(bool)(vars.count("no-vis"));
  int doms;
  double duration;
  double step;
  int order;
  vars.count("dom") ? doms = vars["dom"].as<int>() : doms=2;
  vars.count("dur") ? duration = vars["dur"].as<double>() : duration = 10;
  vars.count("step") ? step = vars["step"].as<double>() : step = .01;
  vars.count("ord") ? order = vars["ord"].as<int>() : order = 20;
  
  //construct the inputs to wave construction

  std::vector<double> x = std::vector<double>();
  std::vector<double> times;
  times.push_back(0.0);

  std::vector<int> orders(doms,order);

  std::vector<elementStateHistory> states;
	  
  std::vector<std::shared_ptr<std::vector<double>>> abscissas(doms);
  std::vector<std::shared_ptr<std::vector<double>>> weights(doms);
  std::vector<std::shared_ptr<matrix<double>>> DMats(doms);
  
  std::transform(orders.begin(),orders.end(),abscissas.begin(),
		 isDG ? legendreTools::generateAbscissas : legendreTools::generateGLAbscissas);
  std::transform(orders.begin(),orders.end(),abscissas.begin(),weights.begin(),
		 isDG ? legendreTools::generateWeights : legendreTools::generateGLWeights);
  std::transform(orders.begin(),orders.end(),abscissas.begin(),DMats.begin(),
		 [](auto ord, auto absc){
		   return legendreTools::generateDMat(ord,absc,legendreTools::generateBaryWeights(ord,absc));} );
  
  
  for(int d=0;d<doms;d++)
    {
      for(int i=0;i<orders[d];i++)
	x.push_back(boundData(-(abscissas[d]->at(i) + 2.0*d )));
      for(int i=0;i<orders[d];i++)
	x.push_back(-boundData(-(abscissas[d]->at(i) + 2.0*d )));
    }

  //Initialize the history with starting scalar functions
  int elstart=0;
  for(int d=0;d<doms;d++)
    {
      states.push_back(elementStateHistory());
      for(int i=0;i<2;i++)
	{
	  states[d].functionStates.push_back(functionStateHistory());
	  states[d].functionStates[i].timeStates.push_back(scalarFunction(orders[d],abscissas[d],weights[d],DMats[d],
									  std::vector<double>(x.begin()+elstart+i*orders[d]
											      ,x.begin()+elstart+(i+1)*orders[d])));
	}
      elstart+=2*orders[d];
    }
  multiStateHistory waveHist(orders,doms,2,states,times);

  //Construct the wave object and evolve it
  if(verb) printf("initializing ode integrator \n");
  if(isDG)
    {
      auto wave = DGTransmittingMultiWave(orders,abscissas,weights,DMats,step,doms,boundData,isReflecting,verb);
      odeEvolve(x,wave,duration,step,waveHist);
    }
  else
    {
      auto wave = collTransmittingMultiWave(orders,abscissas,weights,DMats,step,doms,boundDatadx,isReflecting,verb);
      odeEvolve(x,wave,duration,step,waveHist);
    }
  
  if(verb) printf("Computing legendre modes (summing quadratures)...\n");


  // Dump of full spectral data in form matching hierarchy of the history object
  if(dumpData)
    {
      printf("--Data dump of scalar wave history--\n");
      for(int d=0;d<doms;d++)
	{
	  printf(" domain %d\n",d);
	  for(int f=0;f<2;f++)
	    {
	      printf("  function %d\n",f);
	      for(int t=0;t<states[0].functionStates[0].timeStates.size();t++)
		{
		  printf("   t=%f\n",times[t]);
		  printf("    ");
		  for(int c=0;c<order-1;c++)
		      printf("%f, ",states[d].functionStates[f].timeStates[t].collocationData[c]);
		  printf("%f\n",states[d].functionStates[f].timeStates[t].collocationData[order-1]);
		}
	    }
	}
    }
  
  if(!vis)
      return 0;
  
  // plot the movie of the wavefunction
  std::vector<std::vector<scalarFunction>> plotAccumulator;
  int stepSize = (duration/(step * 1000)) + 1;
  for(int i = 0;i<states[0].functionStates[0].timeStates.size()-stepSize + 1;i+=stepSize)
    {
      plotAccumulator.push_back(std::vector<scalarFunction>());
      for(int j=0;j<states.size();j++)
	plotAccumulator[i/(stepSize)].push_back(states[j].functionStates[0].timeStates[i]);
      multiPlotWaveandDeriv(plotAccumulator[i/stepSize]);
    }

  //plot the top 3 wavemodes in each domain as a function of time

  
  multiPlotTopNModes(plotAccumulator,doms,3,duration,duration/(step*10000) + 1);

  //wait a moment
  sleep(3);

  //plot the bottom 3 wavemodes in each domain as a function of time
  multiPlotBottomNModes(plotAccumulator,doms,3,duration,duration/(step*10000) + 1);
  return 0;
}
  

