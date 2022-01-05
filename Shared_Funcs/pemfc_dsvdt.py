""" Import needed modules """
"-----------------------------------------------------------------------------"
import sys
import numpy as np
import pandas as pd
import cantera as ct
from math import tanh
from scipy.integrate import solve_ivp
from Shared_Funcs.pemfc_transport_funcs import *

""" Define CL dsvdt for core-shell model """
"-----------------------------------------------------------------------------"
def dsvdt_cl(t, sv, dsvdt, gdl_BC):    
    
    eps_w_cutoff = 1e-10
    eps_w_dropoff = 1e-5
        
    """ Set up conditions at GDL/CL BC """    
    # Load in BC state and flux from GDL:--------------------------------------
    flux_g_up = gdl_BC['flux_g_up']
    flux_w_up = gdl_BC['flux_w_up']
    
    rho_gas_k = sv[ca.ptr['rho_gas_k']]
    TDY1 = cl.d['T'], sum(rho_gas_k), rho_gas_k
    eps_w1 = sv[ca.ptr['eps_w_cl']]
    
    i_io_up = 0  # no protons flow into the GDL
        
    """ Begin loop - with if statements for CL/Elyte BC """
    for i in range(cl.d['Ny']):                
        # Set the states for outer reactions:----------------------------------
        ca.outter_rxn_state(cl,sv,i)
        
        # Smooth step change for rxns dependent on h2o surface area:
        eps_w = max(sv[ca.ptr['eps_w_cl'] +i*cl.d['nxt_y']], eps_w_cutoff)
        mult = tanh(eps_w / eps_w_dropoff)
        
        # Gas-phase production/consumption terms:            
        sdot_ng_gk = ca.naf_gs[i].get_net_production_rates(ca.gas)*cl.d['SApv_ng'][i]
        sdot_wg_gk_p = ca.h2o_gs.get_creation_rates(ca.gas)*cl.d['SApv_wg'][i]*mult
        sdot_wg_gk_d = ca.h2o_gs.get_destruction_rates(ca.gas)*cl.d['SApv_wg'][i]
        
        rho_dot_ng_gk = sdot_ng_gk*cl.d['1/eps_g'][i]*ca.gas.molecular_weights
        rho_dot_wg_gk = (sdot_wg_gk_p - sdot_wg_gk_d)*cl.d['1/eps_g'][i]\
                      *ca.gas.molecular_weights
                      
        # Nafion phase production/consumption terms:
        sdot_ng_nk = ca.naf_gs[i].get_net_production_rates(ca.naf_b[i])*cl.d['SApv_ng'][i]
        sdot_nw_nk_p = ca.naf_ws[i].get_creation_rates(ca.naf_b[i])*cl.d['SApv_nw'][i]*mult
        sdot_nw_nk_d = ca.naf_ws[i].get_destruction_rates(ca.naf_b[i])*cl.d['SApv_nw'][i]
        
        rho_dot_ng_nk = sdot_ng_nk*cl.d['1/eps_n'][i]*ca.naf_b[i].molecular_weights
        rho_dot_nw_nk = (sdot_nw_nk_p - sdot_nw_nk_d)*cl.d['1/eps_n'][i]\
                      *ca.naf_b[i].molecular_weights
                   
        # H2O(Liq) production/consumption terms:
        sdot_nw_w_p = ca.naf_ws[i].get_creation_rates(ca.h2o_b)*cl.d['SApv_nw'][i]
        sdot_nw_w_d = ca.naf_ws[i].get_destruction_rates(ca.h2o_b)*cl.d['SApv_nw'][i]*mult
        sdot_wg_w_p = ca.h2o_gs.get_creation_rates(ca.h2o_b)*cl.d['SApv_wg'][i]
        sdot_wg_w_d = ca.h2o_gs.get_destruction_rates(ca.h2o_b)*cl.d['SApv_wg'][i]*mult
        
        eps_dot_nw_w = (sdot_nw_w_p - sdot_nw_w_d)*ca.h2o_b.volume_mole
        eps_dot_wg_w = (sdot_wg_w_p - sdot_wg_w_d)*ca.h2o_b.volume_mole
                  
        # Set states for inner reactions (ORR):--------------------------------
        ca.inner_rxn_state(cl,sv,i)
        
        sdot_nc_e = ca.pt_s[i].get_net_production_rates(ca.carb) 
        sdot_nc_pt = ca.pt_s[i].get_net_production_rates(ca.pt_s[i])
        sdot_nc_nk = ca.pt_s[i].get_net_production_rates(ca.naf_b[i])
        
        i_Far = sdot_nc_e *ct.faraday *cl.d['SApv_pt'][i]
        
        theta_dot_pt = sdot_nc_pt *cl.d['1/gamma']
        
        rho_dot_nc_nk = sdot_nc_nk *cl.d['SApv_pt'][i] *cl.d['1/eps_n'][i]\
                      *ca.naf_b[i].molecular_weights
        
        # Double layer potential at each Y node:-------------------------------

        # BC for CL and electrolyte interface
        if i == cl.d['Ny'] -1: 
            i_io_dwn = cl.d['i_ext']
            
        # Internal CL nodes
        else:
            sig_io = np.mean(cl.d['sig_io'][i:i+2])
            epstau_n = np.mean(cl.d['eps/tau2_n'][i:i+2])
            
            phi_elyte1 = sv[ca.ptr['phi_dl'] +i*cl.d['nxt_y']]
            phi_elyte2 = sv[ca.ptr['phi_dl'] +(i+1)*cl.d['nxt_y']]
            
            i_io_dwn = sig_io *epstau_n *(phi_elyte1 - phi_elyte2) *cl.d['1/dy']        

        i_dl = (i_io_up - i_io_dwn) *cl.d['1/dy'] - i_Far 
        dsvdt[ca.ptr['phi_dl'] +i*cl.d['nxt_y']] = i_dl *cl.d['1/CA_dl'][i]

        i_io_up = i_io_dwn
        
        # V_w and Gas phase species at each Y node:----------------------------
        
        # BC for CL and electrolyte interface
        if i == cl.d['Ny'] -1: 
            flux_g_dwn = np.zeros(ca.gas.n_species)
            flux_w_dwn = 0 
            
        # Internal CL nodes
        else:
            rho_gas_k = sv[ca.ptr['rho_gas_k'] +(i+1)*cl.d['nxt_y']]
            TDY2 = cl.d['T'], sum(rho_gas_k), rho_gas_k
            
            eps_w2 = sv[ca.ptr['eps_w_cl'] +(i+1)*cl.d['nxt_y']]
            TDY_vec, eps_vec = [TDY1, TDY2], [eps_w1, eps_w2]
            
            flux_g_dwn = fickian_adf(TDY_vec, ca, cl, cl.d, cl.d, i)
            flux_w_dwn = darcys_law(TDY_vec, eps_vec, ca, cl, cl.d, cl.d, i)
          
        # Water volume fraction time evolution
        w_i = ca.ptr['eps_w_cl'] +i*cl.d['nxt_y']     
        dsvdt[w_i] = (flux_w_up - flux_w_dwn) *cl.d['1/dy'] *ca.h2o_b.volume_mole\
                    + eps_dot_nw_w +eps_dot_wg_w 
                      
        # Gas-phase densities time evolution  
        g_i = ca.ptr['rho_gas_k'] +i*cl.d['nxt_y']
        dsvdt[g_i] = (flux_g_up - flux_g_dwn) *cl.d['1/dy'] *cl.d['1/eps_g'][i]\
                   + rho_dot_ng_gk +rho_dot_wg_gk +sv[g_i] *dsvdt[w_i] *cl.d['1/eps_g'][i]

        flux_g_up = flux_g_dwn
        TDY1 = TDY2
        
        flux_w_up = flux_w_dwn
        eps_w1 = eps_w2
        
        # Pt coverages at each Y node:-----------------------------------------
        c_i = ca.ptr['theta_pt_k'] +i*cl.d['nxt_y']
        dsvdt[c_i] = theta_dot_pt    
        
        # Nafion densities at each R node:-------------------------------------
        # The Nafion densities change due to reactions at the outter and inner
        # most shells as well as fluxes between adjacent shells. The direction
        # of storage for the radial terms are done from the outermost shell
        # to the innermost one.
        
        " Start by evaluating the outermost shell "
        rho_k1 = sv[ca.ptr['rho_naf_k'] +i*cl.d['nxt_y']]
        rho_k2 = sv[ca.ptr['rho_naf_k'] +i*cl.d['nxt_y'] +cl.d['nxt_r']]
        rho_flx_inr = radial_fdiff(rho_k1, rho_k2, cl, [i,0])
        
        n_i = ca.ptr['rho_naf_k'] +i*cl.d['nxt_y']
        dsvdt[n_i] = (rho_dot_ng_nk + rho_dot_nw_nk) *cl.d['1/Vf_shl'][i,0]\
                   - rho_flx_inr *cl.d['1/r_j'][i,0]**2 *cl.d['1/t_shl'][i,0]
        
        # Ensure constant proton density                                
        dsvdt[n_i[cl.d['iH_n']]] = 0  

        rho_flx_otr = rho_flx_inr
        rho_k1 = rho_k2

        " Evaluate the inner shell nodes "
        for j in range(1, cl.d['Nr'] -1):
            rho_k2 = sv[ca.ptr['rho_naf_k'] +i*cl.d['nxt_y'] +(j+1)*cl.d['nxt_r']]
            rho_flx_inr = radial_fdiff(rho_k1, rho_k2, cl, [i,j])

            n_i = ca.ptr['rho_naf_k'] +i*cl.d['nxt_y'] +j*cl.d['nxt_r']
            dsvdt[n_i] = (rho_flx_otr - rho_flx_inr)\
                       *cl.d['1/r_j'][i,j]**2 *cl.d['1/t_shl'][i,j]
            
            # Ensure constant proton density
            dsvdt[n_i[cl.d['iH_n']]] = 0 
            
            rho_flx_otr = rho_flx_inr
            rho_k1 = rho_k2                                      
        
        # Innermost Nafion node densities:
        n_i = ca.ptr['rho_naf_k'] +i*cl.d['nxt_y'] +(cl.d['Nr'] -1)*cl.d['nxt_r']
        dsvdt[n_i] = rho_dot_nc_nk *cl.d['1/Vf_shl'][i,-1]\
                   + rho_flx_otr *cl.d['1/r_j'][i,-1]**2 *cl.d['1/t_shl'][i,-1]
                   
        # Ensure constant proton density
        dsvdt[n_i[cl.d['iH_n']]] = 0 
        
    return dsvdt

""" Define dsvdt for pemfc models - combined GDL and CL """
"-----------------------------------------------------------------------------"
def dsvdt_func(t, sv):
    # Initialize output:-------------------------------------------------------       
    dsvdt = np.zeros_like(sv)
    
    # Update porosity values to consider water volume fractions
    gdl.update(ca,sv)
    cl.update(ca,sv)
    
    """ Bondary Condition - GDL and gas gas channel """
    # Densities/Temp of GDL gas species and gas channel (top):-----------------
    ca.gas.TPX = gdl.d['TPX_BC']
    TDY1 = ca.gas.TDY 
    eps_w1 = 0
    d1 = gdl.d

    rho_gdl_k = sv[ca.ptr['rho_gdl_k']]
    TDY2 = gdl.d['T'], sum(rho_gdl_k), rho_gdl_k
    eps_w2 = sv[ca.ptr['eps_w_gdl']]
    d2 = gdl.d
    
    TDY_vec = [TDY1, TDY2]
    eps_vec = [eps_w1, eps_w2]
    
    flux_g_up = fickian_adf(TDY_vec, ca, gdl, d1, d2, 'channel')
    flux_w_up = darcys_law(TDY_vec, eps_vec, ca, gdl, d1, d2, 'channel')
    
    TDY1 = TDY2
    eps_w1 = eps_w2
    d1 = d2
    
    for i in range(gdl.d['Ny']):     
        
        # At the GDL/CL boundary (bottom):
        if i == gdl.d['Ny']-1:
            rho_gas_k = sv[ca.ptr['rho_gas_k']]
            TDY2 = cl.d['T'], sum(rho_gas_k), rho_gas_k
            
            eps_w2 = sv[ca.ptr['eps_w_cl']]
            d, d2 = gdl_cl, cl.d
            
        # Internal GDL nodes:
        else:
            rho_gdl_k = sv[ca.ptr['rho_gdl_k'] +(i+1)*gdl.d['nxt_y']]
            TDY2 = gdl.d['T'], sum(rho_gdl_k), rho_gdl_k
            
            eps_w2 = sv[ca.ptr['eps_w_gdl'] +(i+1)*gdl.d['nxt_y']]
            d, d2 = gdl, gdl.d
            
        TDY_vec = [TDY1, TDY2]
        eps_vec = [eps_w1, eps_w2]
        
        flux_g_dwn = fickian_adf(TDY_vec, ca, d, d1, d2, i)     
        flux_w_dwn = darcys_law(TDY_vec, eps_vec, ca, d, d1, d2, i)
                    
        # Water volume fraction time evolution
        w_i = ca.ptr['eps_w_gdl'] +i*gdl.d['nxt_y']
        dsvdt[w_i] = (flux_w_up - flux_w_dwn) *gdl.d['1/dy'] *ca.h2o_b.volume_mole     
        
        # Gas-phase densities time evolution
        g_i = ca.ptr['rho_gdl_k'] +i*gdl.d['nxt_y']
        dsvdt[g_i] = (flux_g_up - flux_g_dwn) *gdl.d['1/dy'] *gdl.d['1/eps_g'][i]\
                   + sv[g_i] *dsvdt[w_i] *gdl.d['1/eps_g'][i]
        
        flux_g_up = flux_g_dwn
        TDY1 = TDY2
                          
        flux_w_up = flux_w_dwn  
        eps_w1 = eps_w2   
        
        d1 = d2
    
    # Store BC values to pass into CL function:
    gdl_BC = {}    
    gdl_BC['flux_g_up'] = flux_g_up
    gdl_BC['flux_w_up'] = flux_w_up

    """ Generic loop for interal CL nodes in y-direction """
    dsvdt = dsvdt_cl(t, sv, dsvdt, gdl_BC)

    # print(t)
    # print(dsvdt)

    # if t > 1e-3:    
    #     user_in = input('"Enter" to continue or "Ctrl+c" to cancel.')   
    #     if user_in == KeyboardInterrupt:
    #         sys.exit(0)

    return dsvdt

""" Use integrator to call dsvdt and solve to SS  """
"-----------------------------------------------------------------------------"   
# Funtion to check values in console: 
def print_chk(v1,v2,v3):
    print('t_f:',v1,'i_ext:',v2,'dphi_ss:',v3)
    return None
    
# Create vectors to store outputs:
i_ext = np.hstack([i_OCV, i_ext0, i_ext1, i_ext2])
dphi_ss = np.zeros_like(i_ext)

# Initialize dataframe for raw solution at each current:
df_i = pd.DataFrame()

# Define common index for last CL node's phi_dl:
iPhi_f = int(ca.ptr['phi_dl'] +(cl.d['Ny']-1)*cl.d['nxt_y'])

for i in range(i_ext.size):
    # Update and convert i_ext: A/cm^2 -> A/m^2
    cl.d['i_ext'] = i_ext[i] *100**2
    
    try: 
        if i_ext[i] == 0: t_fin = t_sim[0]
        else: t_fin = t_sim[1]
        
        sol = solve_ivp(lambda t, sv: dsvdt_func(t, sv), [0, t_fin], SV0, 
              method=method, atol=atol, rtol=rtol, jac_sparsity=sparsity)
        
        # Skip loop below if not converged at OCV
        if t_fin - sol.t[-1] > 0.1*t_fin: 
            print('Did not converge in time')
            break
        
        # Calculate extra PEM resistance terms to subtract off:
        sig_io_eff = cl.d['sig_io'][-1] #*cl.d['eps/tau2_n'][-1]
        R_naf = i_ext[i]*(pem.d['R_naf'] +0.5*cl.d['dy'] /sig_io_eff *100**2)
        
        # Store solution and update initial values:
        SV0 = sol.y[:,-1]
        df_temp = pd.DataFrame([np.append(i_ext[i], sol.y[:,-1])])
        df_i = df_i.append(df_temp, ignore_index=True)
        dphi_ss[i] = sol.y[iPhi_f, -1] - dphi_eq_an - R_naf
                        
        print_chk(sol.t[-1],round(cl.d['i_ext']*1e-4,6),round(dphi_ss[i],3))
    
        if i_OCV != 0 or polar == 'no':
            break

    except: 
        print('Failed in dsvdt loop')
        break

# Build dataframe to store processed info
df_p = pd.DataFrame()
df_p['i_ext [A/cm2]'] = i_ext 
df_p['Voltage [V]'] = dphi_ss
