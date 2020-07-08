# ToyMC_Onium

Code to generate a Toy MC from a parametrization of s^star and cos(alpha) and using partonic PDFs from LHAPDF. Currently generates for a given state and in both
energies sqrt(s) = {7, 13}. The distributions can then be used in a fit of their normalizations and relative contribution of parameter beta.

Files and their usage:
  
* qqbarMC.C: code to generate the Toy MC. All variables (particle mass, parameters of the generating distributions, sqrt(s), etc) must be changed inside the code.
The code is called using cycle.sh in order to include the LHAPDF library (current path only works in my computer, will have to be changed to work in others). 
  
  - Creates a file "MC_res.root", which is then renamed MC_res_[fname].root (fname given by user) by cycle.sh
  
* plot_twos.C: code to plot toy MC histograms vs data graphs of the pT/M distribution in different y bins and the y distribution in different pT bins. 
State to be plotted must be defined inside the code. Input is a file of the type MC_res_[state]_beta[beta_value].root

  - Creates a file "MC_vs_Data_[state].root" (state given by the code based on what quarkonium state was plotted) containing the resulting histograms
  and graphs, as well as a TTree with the normalizations applied

* [plot/pull]_codes.C: auxiliary files imported in the main file plot_twos.C. Contain the macros that fill the data graphs with all correction factors applied and
initialize the MC histograms (plot_codes.C) and the macros that receive the data and MC distributions and determine the pulls (pull_codes.C).

* fitNorm.C: takes the result of the plotting and performs a fit, with the normalization and relative contribution of beta as free parameters. This step is done to
test how well the generated distributions match up to data, without having to fix the above parameters by hand. The states to be used in the fit must be changed inside
the code. Input is a set of files of the type MC_vs_Data_[state].root for each state considered.

    - Creates a set of plots giving the pT/M distribution in all y bins for each state, as well as the summarized pulls and relative deviation (all states and energies) 
  and all the relative deviations (for each state and energy)
  
    - Also creates a root file Fit_results.root with all the fit results and fitted plots
    
  * dev_code.C: auxiliary file imported in fitNorm.C. Contains the macro that calculates the relative deviations.
