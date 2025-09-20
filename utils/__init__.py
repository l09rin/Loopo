"""
LOOPO utilities package.
Contains four different scripts:
 * loopOnTraj : actually perform sequential or parallel parsing and analysis of configurations.
* generate_tsteps_cycle.py : to easily generate files with timesteps having linear or logaritmic interval spacings.
* msd.py : calculates the mean-squared displacement on a trajectory in several different formats; it automatically detects the timesteps intervals for averages.
* msd_parallel.py : same as msd.py, thread-parallelized by using multithreading module.
"""
