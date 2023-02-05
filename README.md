# C Interstellar Medium SPH

A baby SPH code for simulating the ISM.

To compile simply:

```make```

To compile a simulation:

```cd sims/sim_name```

```make```

To run a simulation:

```cd sims/sim_name```

```./sim_name.exe```

To clean up a simulation folder or the main library:

```make clean```

Examine the `.config` files for details of each simulation.

Sample analysis is given in `analysis/cisms_analysis.ipynb`. 

Note that for MacOS, you probably have to change the default `gcc` compiler to something else.

IE: `gcc` -> `gcc-12`
