""" Model Description """
"-----------------------------------------------------------------------------"
"""This model is a half cell model of of PEM fuel cell cathode. The model uses 
a core-shell type catalyst layer geometry. The core-shell geometry involves a 
Pt-covered carbon core covered with a thin shell of Nafion."""

""" Model Assumptions """
"-----------------------------------------------------------------------------"
" 1.) Gas properties are constant -> parameters"
" 2.) Temperature is constant -> parameter"
" 3.) Variables are Phi_dl [V] -> only difference important,"
"     Nafion/gas species mass densities [kg/m^3], and surface coverages [-]"
" 4.) The proton density in Nafion is constant due to its structure"
" 5.) The C/Pt phases at every node have the same potential"
" 6.) The ionic conductivity of Nafion is a function of its RH, temperature,"
"     and thickness. It can be modeled via various sub models described in"
"     detail in the pemfc_property_funcs file."

# Note: Interpolations can only be done within the ranges given by the above
#       mentioned paper [1]; therefore, this model only works for properties of 
#       Nafion RH, temperature, and thickness of 20-95 %, 30-60 C, and 4-300 nm
#       respectively. Values outside of these ranges with result in an unknown
#       value of ionic conductivity.

# Note: Another valuable paper for model inputs is "An improved two-dimensional 
#       agglomerate cathode model to study the influence of catalyst layer 
#       structural parameters" by Sun, Peppley, and Karan from 2005 [2]. This
#       includes GDL and CL thickness, porosities, and more.

# Note: Permeability values taken as scaled values from Zenyuk, Das, Weber 
#       paper that reported saturated values for GDL and CL at reported
#       porosities. Title "Understanding Impacts of Catalyst Layer Thickness
#       on Fuel Cell Performance via Mathematical Modeling" (2016) [3]. 

# Note: Oxygen diffusion through Nafion was taken from Sethuraman et al. paper 
#       and scaled by water volume fraction as a function of t_naf. Titled 
#       "Measuring Oxygen, Carbon Monoxide, Hydrogen Sulfide Diffusion 
#       Coefficients and Solubility in Nafion Membranes" (2009) [4].

# Note: Knowledge of O2 diffusion scaling with water volume fraction used to
#       scale values. Interpolations from DeCaluwe et al. taken to evaluate
#       water fraction as function of t_naf. Titled "Structure-property 
#       relationships at Nafion thin-film interfaces:..." (2018) [5]

# Note: Nafion conducitivity treatment as bulk material was added as an option.
#       This method assumes the thin shell has the same conductivity as bulk
#       material and that there are no confinement or microstructure effects
#       when it gets too thin. The relationship is taken from Yadav et al. in
#       their 2012 publication of "Analysis of EIS Technique and Nafion 117 
#       Conductivity as a Function of Temperature and Relative Humidity" [6]

# Note: Low Pt loading data and modeling results are present with operating 
#       conditions in "Modeling and Experimental Validation of Pt Loading and
#       Electrode Composition Efects in PEM Fuel Cells" by L. Hao et. al. [7]
#       This paper was published in 2015 and makes it recent enough to validate
#       against.

# Note: In order to validate the model, results for the full cell are needed.
#       The anode was added by determining the steady-state potential via
#       Gibb's free energy correlations. The PEM was simply modeled as a V=IR
#       relationship where the resistance (Ohm*cm^2) was taken from [8] (2002).

""" Import needed modules """
"-----------------------------------------------------------------------------"
import os, time
import numpy as np
import cantera as ct

""" User Input Parameters """
"-----------------------------------------------------------------------------"
" External currents "
i_OCV  = 0.                       # 0 [A/cm^2] -> or single run if != 0 
i_ext0 = np.linspace(1e-3,0.1,10) # currents in kinetic region [A/cm^2]
i_ext1 = np.linspace(2e-1,1.5,10) # currents in ohmic region [A/cm^2]
i_ext2 = np.linspace(1.80,2.1,10) # currents near limitting [A/cm^2]

" Transport and material properties "
C_dl = 1.5e-9                   # capacitance of double layer [F/m^2]
sig_meth = 'mix'                # 'lam', 'bulk', 'mix', or 'lit' for sig_io
D_O2_meth = 'mix'               # 'lam', 'bulk', 'mix', or 'sun' for D_O2
R_naf = 0.028                   # ASR of Nafion membrane [Ohm*cm^2]
st_h2o = 6.25e-2                # water surface tension [N/m]
r_p_cl = 25e-9                  # pore radius for cl capillary pressure [m]
r_p_gdl = 15e-6                 # pore radius for gdl capillary pressure [m]     
ca_cl = 110                     # cl water contact angle [degrees]
ca_gdl = 135                    # gdl water contact angle [degrees]
mu_w = 4.97e-4                  # liquid water viscosity [Pa*s]

" Pt loading and geometric values "
w_Pt = 0.2                      # Pt loading per geometric area [mg/cm^2]
rho_Pt = 21.45e3                # Pt density for area property calcs [kg/m^3]
r_c = 50e-9                     # radius of single carbon particle [m] 
r_Pt = 1e-9                     # radius of Pt 1/2 sphere sitting on C [m]
t_gdl = 250e-6                  # thickness of cathode GDL [m]
t_cl = 15e-6                    # thickness of cathode CL [m]
t_ref = t_cl                    # CL reference thickness [m]
eps_g_gdl = 0.5                 # porosity of cathode GDL [-]
eps_g_cl = 0.1                  # porosity of cathode CL [-]
theta = 45                      # O2 transport angle (<90) [degrees]
t_naf = np.ones(5)*12e-9        # thickness of Nafion shells [m]
pt_dist = 'equal'               # CL distribution of pt: 'lin', 'exp', 'equal'

" Model inputs "
method = 'BDF'                  # method for solve_ivp [eg:BDF,LSODA,...]
t_sim  = [1e6,1e2]              # time span of integration [s]
Ny_gdl = 10                     # number of depth nodes for GDL
Ny_cl = 5                       # number of depth nodes for CL
Nr_cl = 3                       # number of radial nodes for CL Nafion shells

" Modify tolerance for convergence "
atol = 1e-9                     # absolute tolerance passed to solver
rtol = 1e-4                     # relative tolerance passed to solver

" Plot toggles - (0:off, 1:on, 2:extra options) "
post_only = 2       # only run post-processing (1:current, 2:past)
debug  = 0          # plot solution vector variables vs time (1:low, 2:high)
grads  = 1          # plot solution vector depth gradients (1:low, 2:high)
radial = 0          # plot solution vecotr radial gradients (1:low, 2:high)
polar  = 1          # plot full cell polarization curve (1:w/o pwr, 2:w/ pwr)
over_p = 0          # plot overpotential curve

" Verification settings "
data = 0            # include data from Owejan et. al. if (0: off, 1: on)
i_ver = 1           # verify current between GDL and CL (0: off, 1: on)
i_find = 1.83       # current to use in i_ver calculation [A/cm^2]

" Save/load options "
load_folder = 'mix_2_air'            # folder name to load past dfs 
save_folder = 'folder_name'          # folder name for saving files
save = 0                             # toggle saving on/off

" Initialization and BC inputs "
phi_ca_init = 0.7                    # initial cathode potential [V]
T_ca, P_ca = 333, 1.5*ct.one_atm     # cathode T [K] and P [Pa] 
T_an, P_an = 333, 1.5*ct.one_atm     # anode T [K] and P [Pa] 

X_BC = 'N2: 0.79,O2: 0.21,H2O: 0.15' # cathode gas channel composition  [-]
T_BC, P_BC = 333, 1.5*ct.one_atm     # cathode gas channel T [K] and P [Pa]

" Reaction properties "
n_elec_an = 2                        # sum(nu_k*z_k) for anode surf rxn


""" End of user inputs - do not edit anything below this line """
"-----------------------------------------------------------------------------"
###############################################################################
###############################################################################

""" Process inputs from this file and run model """
"-----------------------------------------------------------------------------"
if __name__ == '__main__':  
    exec(open("Shared_Funcs/pemfc_pre.py").read())
    
    if post_only == 0:
        time1 = time.time()
        exec(open("Shared_Funcs/pemfc_dsvdt.py").read())
        time2 = time.time()
        print('\nRun Time: ', np.round((time2-time1)/60,2),' min')
    
    exec(open("Shared_Funcs/pemfc_post.py").read())    
    
# chi2 for air runs:
# Current mix best fit is chi2 = 1.90
# Current bulk best fit is chi2 = 2.71
# Current lit relation (w/ RH dependence on gas) best fit is chi2 = 3.92

# chi2 for just w_Pt = 0.025: 
    # 2.52, 2.33, 5.08
    
# chi2 for o2 runs:
# Current mix best fit is chi2 = 2.61
# Current bulk best fit is chi2 = 2.62
# Current lit relation (w/ RH dependence on gas) best fit is chi2 = 5.06
    
# Explanation of above:
# The mix and bulk are as described by the paper and pemfc_property_funcs file.
# As for the literature values, the w/ dependence oon rho_naf_h2o uses the 
# RH_eq function, which still takes into consideration reduced water uptake. On
# the other hand, the w/ only dependence on gas RH actually uses the local gas 
# RH to calculate the conductivity for the CL Nafion. 
    
def chk(sv):
    
    io2_n = ca.naf_b[0].species_index('O2(Naf)')
    io2_g = ca.gas.species_index('O2')
    
    rho_naf_k = sv[ca.ptr['rho_naf_k']]
    rho_gas_k = sv[ca.ptr['rho_gas_k']]
    
    ca.naf_b[0].Y = rho_naf_k
    ca.gas.TDY = cl.d['T'], sum(rho_gas_k), rho_gas_k
    
    c_k_n = ca.naf_b[0].density_mole*ca.naf_b[0].X
    c_k_g = ca.gas.density_mole*ca.gas.X
    
    c_o2_n = c_k_n[io2_n]
    c_o2_g = c_k_g[io2_g]
    
    return print('C_O2(Naf): ', c_o2_n, 'C_O2(gas): ', c_o2_g)
    
    