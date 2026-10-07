This repository includes MATLAB code and data to accompany the manuscript 'Groundwater feedbacks on ice sheets and subglacial hydrology' by GJ Cairns, GP Benham and IJ Hewitt. Code was written by GJ Cairns, with the variable meshgrid in 'variable_grid_300.mat' produced from code by IJ Hewitt.

The script 'make_steady_state_figures.m' can be run to generate the data in 'steady states/', and to create Figures 2, 3, 4 and 5 of the manuscript. The script 'make_steady_state_figures.m' can be run to generate the data in 'retreat_advance/' and create Figures 6, A1 and A2. Other scripts are auxiliary functions used in the above:
* 'default_values.m' generates default values for parameters, scalings and functions used repeatedly in the solutions,
* 'find_steady_state.m' generates steady state solutions, including the data in 'steady_states/' and the initial/final states of advance and retreat,
* 'make_rough_steady_state.m' uses an asymptotic solution to generate an initial guess for the steady state solver,
* 'OBJ_FUNC_ice_hyd.m' is the objective function to solve the coupled ice and hydrology equations at a certain timestep,
* 'retreat_or_advance.m' generates the solutions involving grounding line advance and retreat stored in 'retreat_advance/',
* 'variable_grid_300.mat' contains the meshgrid on which all solutions are found.

The governing equations are solved for *H*, *U*, *x_g* and *N* at a given time on a stretched grid. A 301-point meshgrid is used, with variable spacing so that meshpoints are denser near the grounding line. *H* and *N* are found at whole points, and *U* at half meshpoints.

Please note that these scripts make use of the *cmocean* colormap package (Thyng et al., 2016) in order to produce figures.

References:
Thyng, K.M., C.A. Greene, R.D. Hetland, H.M. Zimmerle, and S.F. DiMarco. 2016. True colors of oceanography: Guidelines for effective and accurate colormap selection. Oceanography 29(3):9–13. http://dx.doi.org/10.5670/oceanog.2016.66
