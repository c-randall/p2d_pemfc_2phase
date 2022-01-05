""" Import needed modules """
"-----------------------------------------------------------------------------"
from scipy.integrate import solve_ivp
from Shared_Funcs.pemfc_transport_funcs import *
import cantera as ct
import numpy as np
import sys

""" Control options for derivative functions """
"-----------------------------------------------------------------------------"
# Toggles to turn on/off in/outer rxns, gas transports:------------------------
# These toggles are only used for debugging purposes. By leaving on specific
# physics in the model, issues can be narrowed down.
pt_rxn = 1
ns_rxn = 1
gas_tog = 1
gdl_tog = 1

# Only perform surface site tracking with a 2-step mechanism:------------------
if rxn_mech == '1s': pt_rxn = 0
elif rxn_mech == '2s': pt_rxn = 1

""" Define CL dsvdt for core-shell model """
"-----------------------------------------------------------------------------"
def dsvdt_cl_cs(t, sv, dsvdt, gdl_BC):
    """ Set up conditions at GDL/CL BC """
    # Initialize indecies for looping:-----------------------------------------
    cl_ymv = 0  # CL y direction mover (y: GDL -> Elyte)
    
    # Load in BC state and flux from GDL:--------------------------------------
    TDY1 = gdl_BC['TDY1']
    flux_up = gdl_BC['flux_up']
    i_io_up = 0  # no protons flow into the GDL
    
    """ Begin loop - with if statements for CL/Elyte BC """
    for i in range(cl['Ny']):
        
        # Set reaction terms:--------------------------------------------------
        # O2 absorption rxn:
        rho_gas_k = sv[ca.ptr['rho_gas_k'] +cl_ymv]
        rho_naf_k = sv[ca.ptr['rho_naf_k'] +cl_ymv]
        
        ca.gas.TDY = cl['T'], sum(rho_gas_k), rho_gas_k
        ca.naf_b.TDY = cl['T'], sum(rho_naf_k), rho_naf_k
        
        rho_dot_go = ca.naf_s.get_net_production_rates(ca.gas) *cl['SApv_naf']\
                   *cl['1/eps_g'] *ca.gas.molecular_weights *gas_tog
        rho_dot_no = ca.naf_s.get_net_production_rates(ca.naf_b) *cl['SApv_naf']\
                   *cl['1/eps_n'] *ca.naf_b.molecular_weights
                  
        # ORR at Pt surface:
        ca.carb.electric_potential = 0
        ca.pt_s.electric_potential = 0

        ca.naf_b.electric_potential = -sv[ca.ptr['phi_dl'] +cl_ymv]
        ca.naf_s.electric_potential = -sv[ca.ptr['phi_dl'] +cl_ymv]
        
        rho_k_inr = sv[ca.ptr['rho_naf_k'] +cl_ymv +(cl['Nr']-1)*cl['nxt_r']]
        ca.pt_s.coverages = sv[ca.ptr['theta_pt_k'] +cl_ymv]
        ca.naf_b.TDY = cl['T'], sum(rho_k_inr), rho_k_inr
        
        rho_dot_gi = ca.pt_s.get_net_production_rates(ca.gas) *cl['SApv_pt']\
                   *cl['1/eps_g'] *ca.gas.molecular_weights *gas_tog
        rho_dot_ni = ca.pt_s.get_net_production_rates(ca.naf_b) *cl['SApv_pt']\
                   *cl['1/eps_n'] *ca.naf_b.molecular_weights
        
        # Gas phase species at each Y node:------------------------------------
        if i == cl['Ny'] -1: # BC for CL and electrolyte interface
            flux_dwn = np.zeros(ca.gas.n_species)
        else:
            rho_gas_k = sv[ca.ptr['rho_gas_k'] +cl_ymv +cl['nxt_y']]
            TDY2 = cl['T'], sum(rho_gas_k), rho_gas_k
            flux_dwn = fickian_adf(TDY1, TDY2, ca.gas, cl, cl, cl, gas_tog)
                  
        # Include O2 rxn, Pt rxn, and fluxes in ODE term:        
        dsvdt[ca.ptr['rho_gas_k'] +cl_ymv] = (flux_up - flux_dwn)*cl['1/eps_g']*cl['1/dy']\
                                        + ns_rxn*rho_dot_go + pt_rxn*rho_dot_gi

        flux_up = flux_dwn
        TDY1 = TDY2
        
        # Nafion densities at each R node:-------------------------------------
        # The Naftion densities change due to reactions at the outter and inner
        # most shells as well as fluxes between adjacent shells. The direction
        # of storage for the radial terms are done from the outermost shell
        # to the innermost one.
        
        " Start by evaluating the outermost shell "
        # This node contains an O2 absorption rxn with the gas phase as well as
        # a maxx flux with the adjacent inner node.
        rho_k1 = sv[ca.ptr['rho_naf_k'] +cl_ymv]
        rho_k2 = sv[ca.ptr['rho_naf_k'] +cl_ymv +cl['nxt_r']]
        rho_flx_inr = radial_fdiff(rho_k1, rho_k2, cl, 0)
        
        # Combine absorption and flux to get overall ODE for Nafion densities:
        dsvdt[ca.ptr['rho_naf_k'] +cl_ymv] = ns_rxn *rho_dot_no *cl['1/Vf_shl'][0]\
                                        - rho_flx_inr *cl['1/r_j'][0]**2 *cl['1/t_shl'][0]
        
        # Ensure constant proton density                                
        dsvdt[ca.ptr['rho_naf_k'][cl['iH_n']] +cl_ymv] = 0  

        rho_flx_otr = rho_flx_inr
        rho_k1 = rho_k2

        " Evaluate the inner shell nodes "
        for j in range(1, cl['Nr'] -1):
            rho_k2 = sv[ca.ptr['rho_naf_k'] +cl_ymv +(j+1)*cl['nxt_r']]
            rho_flx_inr = radial_fdiff(rho_k1, rho_k2, cl, j)

            iMid = ca.ptr['rho_naf_k'] +cl_ymv +j*cl['nxt_r']
            dsvdt[iMid] = (rho_flx_otr - rho_flx_inr) *cl['1/r_j'][j]**2 *cl['1/t_shl'][j]
            
            # Ensure constant proton density
            dsvdt[ca.ptr['rho_naf_k'][cl['iH_n']] +cl_ymv +j*cl['nxt_r']] = 0 
            
            rho_flx_otr = rho_flx_inr
            rho_k1 = rho_k2                                      

        # Pt surface coverages:
        dsvdt[ca.ptr['theta_pt_k'] +cl_ymv] = ca.pt_s.get_net_production_rates(ca.pt_s)\
                                         *cl['1/gamma'] *pt_rxn
        
        # Innermost Nafion node densities:
        iLast = ca.ptr['rho_naf_k'] +cl_ymv +(cl['Nr'] -1)*cl['nxt_r']
        dsvdt[iLast] = pt_rxn *rho_dot_ni *cl['1/Vf_shl'][-1] \
                     + rho_flx_otr *cl['1/r_j'][-1]**2 *cl['1/t_shl'][-1]
        
        # Double layer potential at each Y node:-------------------------------
        # The double layer potential is only stored as a function of CL depth.
        # This means that no local potential gradients are shored in the radial
        # direction throughout the Nafion shells.

        # Find ionic currents and define ODE for phi_dl:
        if i == cl['Ny'] -1: # BC for CL and electrolyte interface
            i_io_dwn = cl['i_ext']
        else:
            i_io_dwn = (sv[ca.ptr['phi_dl'] +cl_ymv] - sv[ca.ptr['phi_dl'] +cl_ymv +cl['nxt_y']])\
                     *cl['sig_naf_io'] *cl['1/dy']

        i_Far = pt_rxn *ca.pt_s.get_net_production_rates(ca.carb) *ct.faraday

        i_dl = (i_io_up - i_io_dwn)*cl['1/dy'] - i_Far*cl['SApv_pt']
        dsvdt[ca.ptr['phi_dl'] +cl_ymv] = i_dl*cl['1/CA_dl']

        i_io_up = i_io_dwn

        # Update Y direction moving index:-------------------------------------
        cl_ymv = cl_ymv +cl['nxt_y']
        
    return dsvdt

""" Define dsvdt for pemfc models - combined GDL and CL """
"-----------------------------------------------------------------------------"
def dsvdt_func(t, sv):
    # Initialize indecies for looping:-----------------------------------------
    gdl_ymv = 0 # GDL y direction mover (y: gas channel -> CL)

    dsvdt = np.zeros_like(sv)
    
    """ Bondary Condition - GDL and gas gas channel """
    # Densities/Temp of GDL gas species and gas channel (top):-----------------
    ca.gas.TPY = gdl['TPY_BC']
    TDY_BC = ca.gas.TDY
    
    # Update porosity to consider water volume fraction
    gdl1 = porosity_func([0], gdl, 0)
    
    # If GDL diffusion is turned on, compare adjacent nodes with ADF flux to 
    # determine the BC composition between the GDL and gas channel.
    rho_gdl_k = sv[ca.ptr['rho_gdl_k']]
    TDY1 = gdl['T'], sum(rho_gdl_k), rho_gdl_k
    flux_up = fickian_adf(TDY_BC, TDY1, ca.gas, gdl, gdl1, gdl1, gdl_tog)
    
    for k in range(gdl['Ny'] -1):
        
        # Update porosity to consider water volume fraction
        gdl2 = porosity_func([0], gdl,0)
        
        # Find the downward flux value
        rho_gdl_k = sv[ca.ptr['rho_gdl_k'] +gdl_ymv +gdl['nxt_y']]
        TDY2 = gdl['T'], sum(rho_gdl_k), rho_gdl_k
        flux_dwn = fickian_adf(TDY1, TDY2, ca.gas, gdl, gdl1, gdl2, gdl_tog)
        
        dsvdt[ca.ptr['rho_gdl_k'] +gdl_ymv] = (flux_up - flux_dwn)*gdl['1/eps_g']*gdl['1/dy']
                                                                                                                                            
        flux_up = flux_dwn
        TDY1 = TDY2
        gdl_ymv = gdl_ymv +gdl['nxt_y']
        
    # Use the composition and state of the last GDL node to calculate the flux
    # into the first CL node.
    rho_gas_k = sv[ca.ptr['rho_gas_k']]
    TDY2 = cl['T'], sum(rho_gas_k), rho_gas_k
    flux_dwn = fickian_adf(TDY1, TDY2, ca.gas, gdl_cl, gdl_cl, gdl_cl, gdl_tog)
    
    dsvdt[ca.ptr['rho_gdl_k'] +gdl_ymv] = (flux_up - flux_dwn)*gdl['1/eps_g']*gdl['1/dy']
    
    flux_up = fickian_adf(TDY1, TDY2, ca.gas, gdl_cl, gdl_cl, gdl_cl, gas_tog)
    TDY1 = TDY2
    
    # Store BC values to pass into CL function:
    gdl_BC = {}
    gdl_BC['TDY1'] = TDY1
    gdl_BC['flux_up'] = flux_up

    """ Generic loop for interal CL nodes in y-direction """
    dsvdt = dsvdt_cl_cs(t, sv, dsvdt, gdl_BC)

#    print(t)
#    print(dsvdt)
#
#    user_in = input('"Enter" to continue or "Ctrl+d" to cancel.')   
#    if user_in == KeyboardInterrupt:
#        sys.exit(0)

    return dsvdt

""" Use integrator to call dsvdt and solve to SS  """
"-----------------------------------------------------------------------------"    
# Create vectors to store outputs:
i_ext = np.hstack([i_OCV, i_ext0, i_ext1, i_ext2])
eta_ss, dphi_ss = np.zeros_like(i_ext), np.zeros_like(i_ext)
sv_save = np.zeros([len(SV0) +1, len(i_ext)])

# Define common index for last CL node's phi_dl:
iPhi_f = int(ca.ptr['phi_dl'] + (Ny_cl-1)*L_cl/Ny_cl)

# Update and convert i_ext: A/cm^2 -> A/m^2
cl['i_ext'] = i_ext[0] *100**2

# try:
sol = solve_ivp(lambda t, sv: dsvdt_func(t, sv), [0, t_sim], SV0, 
                method=method, atol=atol, rtol=rtol, max_step=max_t)
    
#     if t_sim - sol.t[-1] > 0.1*t_sim: # don't store if not developed in time
#         sol.y[:,-1] = np.zeros_like(sol.y[:,-1])
#         i_OCV = None # skip rest of currents if didn't converge for OCV

# except: 
#     i_OCV = None
#     class sol:
#         t = np.zeros(2)
#         y = np.zeros([len(SV0),2])

# # Calculate extra PEM resistance terms to subtract off:
# R_naf_vec = i_ext*(pem['R_naf'] + 0.5*cl['dy'] / cl['sig_naf_io'] *100**2)

# # Store solution and update initial values:
# SV0, sv_save[:,0] = sol.y[:,-1], np.append(i_ext[0], sol.y[:,-1])
# dphi_ss[0] = sol.y[iPhi_f, -1] - dphi_eq_an - R_naf_vec[0]
                
# print('t_f:',sol.t[-1],'i_ext:',round(cl['i_ext']*1e-4,3), 'dPhi:',round(dphi_ss[0],3))

# for i in range(len(i_ext) -1):
#     # Don't run the for loop if i_OCV was not set to 0...
#     if any([all([i == 0, i_OCV != 0]), polar == 'off']): 
#         break
    
#     # Update and convert i_ext: A/cm^2 -> A/m^2
#     cl['i_ext'] = i_ext[i+1] *100**2
    
#     try: 
#         sol = solve_ivp(lambda t, sv: dsvdt_func(t, sv), [0, t_sim], SV0, 
#                         method=method, atol=atol, rtol=rtol, max_step=max_t)
    
#         if t_sim - sol.t[-1] > 0.1*t_sim: # don't store if not developed in time
#             sol.y[:,-1] = np.zeros_like(sol.y[:,-1])
#             break
            
#     except: break
    
#     # Store solution and update initial values:
#     SV0, sv_save[:,i+1] = sol.y[:,-1], np.append(i_ext[i+1], sol.y[:,-1])

#     eta_ss[i+1] = dphi_ss[0] - sol.y[iPhi_f,-1]
#     dphi_ss[i+1] = sol.y[iPhi_f,-1] - dphi_eq_an - R_naf_vec[i+1]

#     print('t_f:',sol.t[-1], 'i_ext:',round(cl['i_ext']*1e-4,3), 'dPhi:',round(dphi_ss[i+1],3))