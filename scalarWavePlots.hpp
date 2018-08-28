#include <stdio.h>
#include <boost/tuple/tuple.hpp>
#include "scalarFunction.hpp"
#include "gnuplot-iostream.h"

#define PLOTRES 300

Gnuplot gp;

void plotSineWave(double simTime)
{
  std::vector<boost::tuple<double,double>> pts;

  double plotMin = 0;
  double plotMax = 10;
  std::printf("current time is %f\n",simTime);
    
  for(double val = plotMin; val<plotMax; val+=(plotMax-plotMin)/PLOTRES) {
    pts.push_back(boost::make_tuple(val,cos(val + simTime)));
  }
  gp << "set term x11 1 noraise\n";
  gp << "set xrange[0:10]\nset yrange [-2:2]\n";
  gp << "plot '-' with lines title 'wave'\n";
  gp.send1d(pts);
  gp << "reread\n";
}


void plotWaveandDeriv(scalarFunction vals)
{
  std::vector<boost::tuple<double,double>> pts;
  std::vector<boost::tuple<double,double>> derivs;
  
  double plotMin = -1;
  double plotMax = 1;

  for(double val= plotMin;val<plotMax;val+=(plotMax-plotMin)/PLOTRES)
    {
      pts.push_back(boost::make_tuple(val,vals.at(val)));
      derivs.push_back(boost::make_tuple(val,vals.dx(val)));
    }
  gp << "set term x11 1 noraise\n";
  gp<< "set xrange[-1:1]\nset yrange[-5:5]\n";
  gp<< "plot '-' with lines title 'simpleWave', '-' with lines title 'deriv'\n";
  gp.send1d(pts);
  gp.send1d(derivs);
  gp << "reread\n";
}

void multiPlotWaveandDeriv(std::vector<scalarFunction> vals)
{
  std::vector<boost::tuple<double,double>> pts;
  std::vector<boost::tuple<double,double>> derivs;

  double plotMin = -1;
  double plotMax = -1 + 2.0*vals.size();

  for(double val= plotMin;val<plotMax;val+=(plotMax-plotMin)/PLOTRES)
    {
      pts.push_back(boost::make_tuple(val,vals[(int)((val+1)/2.0)].at(val - 2.0 * (int)((val + 1)/2.0)) ));
      derivs.push_back(boost::make_tuple(val,vals[(int)((val+1)/2.0)].dx(val - 2.0 * (int)((val + 1)/2.0)) ));
    }
  gp << "set term x11 1 noraise\n";
  gp<< "set xrange[-1:"<<plotMax<<"]\nset yrange[-5:5]\n";
  gp<< "plot '-' with lines title 'simpleWave', '-' with lines title 'deriv'\n";
  gp.send1d(pts);
  gp.send1d(derivs);
  gp << "reread\n";
}

void multiPlotTopNModes(std::vector<std::vector<scalarFunction>> vals,int doms, int n,int maxtime, int timesteps)
{
  std::vector<std::vector<boost::tuple<double,double>>> modes;
      
  double plotMin = 0;
  double plotMax = maxtime;

  for(int d=0;d<doms;d++)
    {
      for(int i=0;i<n;i++)
	{
	  modes.push_back(std::vector<boost::tuple<double,double>>());
	  for(int j=0;j<vals.size();j+=timesteps)
	    {
	      vals[d].at(j).quadSum();
	      modes[d*n+i].push_back(boost::make_tuple(j,vals[d].at(j).spectralData->at(n-1-i)));
	    }
	}
    }
  gp<< "set xrange[0:"<<maxtime<<"]\nset yrange[-10:10]\n";
  gp<<"plot ";
  for(int d=0;d<doms;d++)
    {
      for(int i=0;i<n-1;i++)
	gp<< "'-' with lines title 'domain"<< d <<", wavemode "<<i<<"',";
      if(d!=doms-1)
	gp<<"'-' with lines title 'domain"<< d  <<", wavemode "<<n-1<<"',";
	
    }
  gp<<"'-' with lines title 'domain"<< doms-1  <<", wavemode "<<n-1<<"'\n";
  for(int i=0;i<n*doms;i++)
    gp.send1d(modes[i]);
  gp << "reread\n";
}


void plotTopNModes(std::vector<scalarFunction> vals,int n,int maxtime, int timesteps)
{
  std::vector<std::vector<boost::tuple<double,double>>> modes;
      
  double plotMin = 0;
  double plotMax = maxtime;

  for(int i = 0; i<n;i++)
    {
      modes.push_back(std::vector<boost::tuple<double,double>>());
      for(int j=0;j<vals.size();j+=timesteps)
	{
	  vals.at(j).quadSum();
	  modes[i].push_back(boost::make_tuple(j,vals.at(j).spectralData->at(n-1-i)));
	}
    }
  gp<< "set xrange[0:"<<maxtime<<"]\nset yrange[-10:10]\n";
  gp<<"plot ";
  for(int i=0;i<n-1;i++)
    gp<< "'-' with lines title 'wavemode "<<i<<"',";
  gp<<"'-' with lines title 'wavemode "<<n-1<<"'\n";
  for(int i=0;i<n;i++)
    gp.send1d(modes[i]);
  gp << "reread\n";
}

void plotWave(scalarFunction vals)
{
  std::vector<boost::tuple<double,double>> pts;

  double plotMin = -1;
  double plotMax = 1;

  for(double val= plotMin;val<plotMax;val+=(plotMax-plotMin)/PLOTRES)
      pts.push_back(boost::make_tuple(val,vals.at(val)));
  gp<< "set xrange[-1:1]\nset yrange[-5:5]\n";
  gp<< "plot '-' with lines title 'simpleWave'\n";
  gp.send1d(pts);
  gp << "reread\n";
}

void plotWaveDx(scalarFunction vals)
{
  std::vector<boost::tuple<double,double>> pts;

  double plotMin = -1;
  double plotMax = 1;

  for(double val= plotMin;val<plotMax;val+=(plotMax-plotMin)/PLOTRES)
      pts.push_back(boost::make_tuple(val,vals.dx(val)));
  gp<< "set xrange[-1:1]\nset yrange[-5:5]\n";
  gp<< "plot '-' with lines title 'simpleWave'\n";
  gp.send1d(pts);
  gp << "reread\n";
}
