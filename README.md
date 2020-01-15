# VSPIN

VSPIN is an integrator for solving the spin evolution of planets (or other satellites). Its main distinguishing feature is that it takes into account effects of companion planets in or near a mean motion resonance. This code is part of the work my PhD thesis work during my time at UCLA. 

Results based on this code (or earlier versions of it) can be found so far in two papers published in the Monthly Notices of the Royal Astronomical Society: Vinson & Hansen (2017) and Vinson, Tamayo, & Hansen (2019)


### Analytic Version (analytic_orbit)

This version solves for orbital element variations (i.e. mean motion) analytically based on a simple Pendulum model as can be found in Murray and Dermott (1999), Planetary Dynamics. It assumes the presence of one companion satellite near a mean motion resonance. 

### Non-analytic Version (read_orbit)

This version reads in orbital variations which can be taken from a separate N-body integrator. This is better for systems with more than one resonant companions (e.g. TRAPPIST-1 system), which can result in chaotic evolutions.
