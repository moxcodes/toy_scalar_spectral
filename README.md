# toy scalar wave spectral simulation

## usage:

`>./run-wave --help

Options:
--help                show this help message
--id arg              specify initial data type (sin,fastsin,pulse)
--bc arg              specify right boundary condition (transmit,reflect)
--data                dump time-series collocation data to stdout
--dom arg             specify number of domains
--dur arg             duration of simulation
--step arg            size of simulation timestep
--type arg            type of spectral simulation (coll,dg)
--ord arg             spectral order
--no-vis              turn off default visualizations
--verbose             turn on periodic status updates during simulation`

Default is to run a simple transmitting wave using discontinuous Galerkin method
with two domains and transmitting boundary conditions and a modest spectral
order (10 for each domain). The result is then plotted in a brief animation,
followed by two displays of the top several modes, followed (after a pause) by
the bottom several modes for each domain.

Different initial data may be chosen among a sinusoid (sin), a high-frequency
sinusoid (fastsin), and a gaussian pulse (pulse).

Boundary conditions may be chosen between transmitting (transmit) and reflecting
(reflect) boundaries on the right-hand side.

The type of simulation may be chosen between continuous (coll) and discontinuous
Galerkin (dg).

The number of domains, duration of simulation, step size, and order of spectral
approximation may all be specified in command line.
