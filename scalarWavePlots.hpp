#include <stdio.h>
#include <boost/tuple/tuple.hpp>
#include "scalarFunction.hpp"
#include "gnuplot-iostream.h"

#define PLOTRES 300

/// This provides a handful of helper functions for plotting the
/// scalarFunction values extracted from the various wave simulations.  All
/// functions use gnuplot-iostream.h to draw the plots and update them.
namespace scalarPlots{

  Gnuplot gp;///< the global gnuplot variable for piping gnuplot commands to.


  /// plots the values and first derivatives of a scalarFunction from -1 to 1 using
  /// gnuplot-iostream. Assumes single domain. Prepares stream for replotting.
  /// \param vals scalarFunction to be plotted
  void plotWaveandDeriv(scalarFunction &vals)
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

  /// Plots the values and first derivatives of a set of scalar functions
  /// plotted end-to-end using gnuplot-iostream. Assumes several domains in
  /// successive scalarFunctions. Prepares stream for replotting.
  /// \param vals the scalarFunctions to be plotted
  void multiPlotWaveandDeriv(std::vector<scalarFunction> &vals)
  {
    std::vector<boost::tuple<double,double>> pts;
    std::vector<boost::tuple<double,double>> derivs;

    double plotMin = -1;
    double plotMax = -1 + 2.0*vals.size();

    for(double val= plotMin;val<plotMax;val+=(plotMax-plotMin)/PLOTRES)
      {
	pts.push_back(boost::make_tuple(val,vals[(int)((val+1)/2.0)].at((double)(val - 2.0 * (int)((val + 1)/2.0)) )));
	derivs.push_back(boost::make_tuple(val,vals[(int)((val+1)/2.0)].dx((double)(val - 2.0 * (int)((val + 1)/2.0)) )));
      }
    gp << "set term x11 1 noraise\n";
    gp<< "set xrange[-1:"<<plotMax<<"]\nset yrange[-5:5]\n";
    gp<< "plot '-' with lines title 'simpleWave', '-' with lines title 'deriv'\n";
    gp.send1d(pts);
    gp.send1d(derivs);
    gp << "reread\n";
  }

  /// Plot the highest set of modes for a set of scalar functions over the time
  /// of the simulation using gnuplot-iostream. Prepares the stream for
  /// replotting. vals[i][d] should be the function in domain d at timestep i.
  /// \param vals the scalar function data input.
  /// \param doms number of domains
  /// \param n number of modes to plot
  /// \param maxtime duration of the simulation
  /// \param timesteps number of timesteps to skip over in each plot step
  void multiPlotTopNModes(std::vector<std::vector<scalarFunction>> &vals,int doms, int n,double maxtime, int timesteps)
  {
    std::vector<std::vector<boost::tuple<double,double>>> modes;
     
    double plotMin = 0;
    double plotMax = maxtime;

    for(int d=0;d<doms;d++)
      {
	for(int i=0;i<n;i++)
	  {
	    modes.push_back(std::vector<boost::tuple<double,double>>());
	    for(int j=0;j<vals.size()-timesteps+1;j+=timesteps)
	      {
		vals[j].at(d).quadSum();
		modes[d*n+i].push_back(boost::make_tuple(maxtime*(double)j/(double)vals.size(),
							 vals[j].at(d).spectralData->at(vals[j].at(d).spectralData->size()-1-i)));
	      }
	  }
      }
    gp<< "set xrange[0:"<<maxtime<<"]\nset yrange[-10:10]\n";
    gp<<"plot ";
    for(int d=0;d<doms;d++)
      {
	for(int i=0;i<n-1;i++)
	  gp<< "'-' with lines title 'domain"<< d <<", wavemode "<<vals[0].at(0).spectralData->size()-1-i<<"',";
	if(d!=doms-1)
	  gp<<"'-' with lines title 'domain"<< d  <<", wavemode "<<vals[0].at(0).spectralData->size()-n<<"',";
	
      }
    gp<<"'-' with lines title 'domain"<< doms-1  <<", wavemode "<<vals[0].at(0).spectralData->size()-n<<"'\n";
    for(int i=0;i<n*doms;i++)
      gp.send1d(modes[i]);
    gp << "reread\n";
  }

  /// Plot the lowest set of modes for a set of scalar functions over the time
  /// of the simulation using gnuplot-iostream. Prepares the stream for
  /// replotting. vals[i][d] should be the function in domain d at timestep i.
  /// \param vals the scalar function data input.
  /// \param doms number of domains
  /// \param n number of modes to plot
  /// \param maxtime duration of the simulation
  /// \param timesteps number of timesteps to skip over in each plot step
  void multiPlotBottomNModes(std::vector<std::vector<scalarFunction>> &vals,int doms, int n,double maxtime, int timesteps)
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
		vals[j].at(d).quadSum();
		modes[d*n+i].push_back(boost::make_tuple(maxtime*(double)j/(double)vals.size(),
							 vals[j].at(d).spectralData->at(i)));
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

  /// Plot the lowest set of modes for a single scalar function over the time of
  /// the simulation using gnuplot-iostream. Prepares steram for
  /// replotting. vals[i] should be the function at timestep i.
  /// \param vals the scalar function data input.
  /// \param n number of modes to plot
  /// \param maxtime duration of the simulation
  /// \param timesteps number of timesteps to skip over in each plot step
  void plotBottomNModes(std::vector<scalarFunction> &vals,int n,int maxtime, int timesteps)
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
	    modes[i].push_back(boost::make_tuple(j,vals.at(j).spectralData->at(i)));
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

  /// plots the values of a scalarFunction from -1 to 1 using
  /// gnuplot-iostream. Assumes single domain. Prepares stream for replotting.
  /// \param vals scalarFunction to be plotted
  void plotWave(scalarFunction &vals)
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

  /// plots the first derivatives of a scalarFunction from -1 to 1 using
  /// gnuplot-iostream. Assumes single domain. Prepares stream for replotting.w
  /// \param vals scalarFunction to be plotted
  void plotWaveDx(scalarFunction &vals)
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
}
