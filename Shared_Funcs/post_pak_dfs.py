""" 
Create dataframes for saving pemfc output. Organization is done in order to 
assist in debugging efforts and plotting generality.
"""

import numpy as np
import pandas as pd
import cantera as ct

def df_debug(sol,ca,gdl,cl):
    # sol: solution object from solve_ivp
    # ca: cathode object with phases and pointer info
    # gdl: GDL dictionary with Ny, Len, nxt_y, y
    # cl: CL dictionary with Ny, Nr, Len, nxt_y, nxt_r, y, 1/r_j
    
    """ Data frame to store variable changes over time: """
    df_t = pd.DataFrame()
    df_t['Time [s]'] = sol.t
    
    # GDL variables
    k, mvr = 'eps_w_gdl', 0
    for i in range(gdl.d['Ny']):
        col = k + ' y' + str(i)
        df_t[col] = sol.y[ca.ptr[k] + mvr,:]
        mvr = mvr + gdl.d['nxt_y']
        
    k = 'rho_gdl_k'
    for i in range(ca.ptr[k].size):
        mvr = 0
        
        for j in range(gdl.d['Ny']):
            col = k + ' y' + str(j) + ' ' + ca.gas.species_names[i]
            df_t[col] = sol.y[ca.ptr[k][i] + mvr,:]
            mvr = mvr + gdl.d['nxt_y']
    
    # CL variables
    k, mvr = 'phi_dl', 0
    for i in range(cl.d['Ny']):
        col = k + ' y' + str(i)
        df_t[col] = sol.y[ca.ptr[k] + mvr,:]
        mvr = mvr + cl.d['nxt_y']
        
    k, mvr = 'eps_w_cl', 0
    for i in range(cl.d['Ny']):
        col = k + ' y' + str(i)
        df_t[col] = sol.y[ca.ptr[k] + mvr,:]
        mvr = mvr + cl.d['nxt_y']
        
    k = 'rho_gas_k'
    for i in range(ca.ptr[k].size):
        mvr = 0
        
        for j in range(cl.d['Ny']):
            col = k + ' y' + str(j) + ' ' + ca.gas.species_names[i]
            df_t[col] = sol.y[ca.ptr[k][i] + mvr,:]
            mvr = mvr + cl.d['nxt_y']
            
    k = 'theta_pt_k'
    for i in range(ca.ptr[k].size):
        mvr = 0
        
        for j in range(cl.d['Ny']):
            col = k + ' y' + str(j) + ' ' + ca.pt_s[0].species_names[i]
            df_t[col] = sol.y[ca.ptr[k][i] + mvr,:]
            mvr = mvr + cl.d['nxt_y']
            
    k = 'rho_naf_k'
    for i in range(cl.d['Ny']):         
        coly = k + ' y' + str(i)
        
        for j in range(ca.ptr[k].size):      
            mvr = i*cl.d['nxt_y']
            for n in range(cl.d['Nr']):
                colr = coly + ' r' + str(n) + ' ' + ca.naf_b[0].species_names[j]
                df_temp = pd.DataFrame({colr: sol.y[ca.ptr[k][j] + mvr,:]})
                df_t = pd.concat([df_t,df_temp],axis=1)
                mvr = mvr + cl.d['nxt_r']
        
    return df_t

def df_tool(df_i,df_p,ca,gdl,cl,i_find,yamlfile):
    
    keys = list(ca.ptr.keys())
    
    """ Edit column titles and rearrange df_i """
    df_f = pd.DataFrame()
    df_f['i_ext [A/cm2]'] = df_i[0]
    
    for i in range(len(keys)):
        
        if keys[i] == 'eps_w_gdl':
            for j in range(gdl.d['Ny']):
                pull = j*gdl.d['nxt_y'] + ca.ptr[keys[i]]+1
                df_f[keys[i]+' y'+str(j)] = df_i[pull]
                
        elif keys[i] == 'rho_gdl_k':
            for k,sp in enumerate(ca.gas.species_names):
                for j in range(gdl.d['Ny']):
                    pull = j*gdl.d['nxt_y'] + ca.ptr[keys[i]][k]+1
                    df_f[keys[i]+' y'+str(j)+' '+sp] = df_i[pull]
                    
        elif keys[i] == 'phi_dl':
            for j in range(cl.d['Ny']):
                pull = j*cl.d['nxt_y'] + ca.ptr[keys[i]]+1
                df_f[keys[i]+' y'+str(j)] = df_i[pull]
        
        elif keys[i] == 'eps_w_cl':
            for j in range(cl.d['Ny']):
                pull = j*cl.d['nxt_y'] + ca.ptr[keys[i]]+1
                df_f[keys[i]+' y'+str(j)] = df_i[pull]
                
        elif keys[i] == 'rho_gas_k':
            for k,sp in enumerate(ca.gas.species_names):
                for j in range(cl.d['Ny']):
                    pull = j*cl.d['nxt_y'] + ca.ptr[keys[i]][k]+1
                    df_f[keys[i]+' y'+str(j)+' '+sp] = df_i[pull]
                    
        elif keys[i] == 'theta_pt_k':
            for k,sp in enumerate(ca.pt_s[0].species_names):
                for j in range(cl.d['Ny']):
                    pull = j*cl.d['nxt_y'] + ca.ptr[keys[i]][k]+1
                    df_f[keys[i]+' y'+str(j)+' '+sp] = df_i[pull]
                    
        elif keys[i] == 'rho_naf_k':
            for j in range(cl.d['Ny']):
                for k,sp in enumerate(ca.naf_b[0].species_names):         
                    for r in range(cl.d['Nr']):      
                        pull = j*cl.d['nxt_y'] + r*cl.d['nxt_r'] + ca.ptr[keys[i]][k]+1
                        df_temp = pd.DataFrame({keys[i]+' y'+str(j)+' r'+str(r)+' '+sp: df_i[pull]})
                        df_f = pd.concat([df_f,df_temp],axis=1)             
    
    """ Add power and overpotential to df_p """
    df_p['Power [W/cm2]'] = df_p['i_ext [A/cm2]']*df_p['Voltage [V]']
    df_p['Eta [V]'] = df_p['Voltage [V]'][0] - df_p['Voltage [V]']
    
    """ Determine row index with i_find """
    df_ind = np.argmin(abs(df_f['i_ext [A/cm2]'] - i_find))
    df_row = pd.DataFrame(df_f.loc[df_ind]).T
    df_row = df_row.reset_index(drop=True)
    
    """ Data frame to store steady-state depth gradients: """
    df_y = pd.DataFrame()
    
    y1_gdl, y2_gdl = 0.5*gdl.d['dy'], gdl.d['y']-0.5*gdl.d['dy']
    y1_cl, y2_cl = gdl.d['y']+0.5*cl.d['dy'], gdl.d['y']+cl.d['y']-0.5*cl.d['dy']
    
    df_y['Depth [um]'] = np.hstack([np.linspace(y1_gdl,y2_gdl,gdl.d['Ny']),
                                    np.linspace(y1_cl,y2_cl,cl.d['Ny'])]) *1e6
    
    # Filler for CL-only variables
    fill = np.zeros([gdl.d['Ny'],1])
    
    # Water volume fractions
    eps_cols = [col for col in df_row.columns if 'eps_w' in col]        
    df_y['eps_w [-]'] = df_row[eps_cols].T.to_numpy()
    
    # Gas phase species' densities
    for sp in ca.gas.species_names:
        sp_cols = [col for col in df_row.columns if sp in col]
        gas_cols = [col for col in sp_cols if not '(Naf)' in col]        
        df_y['rho_k '+sp+'(gas)'] = df_row[gas_cols].T.to_numpy()
        
    # Double layer potential
    phi_cols = [col for col in df_row.columns if 'phi_dl' in col]        
    df_y['phi_dl'] = np.vstack([fill,df_row[phi_cols].T.to_numpy()])
    
    # Pt surface coverages
    for sp in ca.pt_s[0].species_names:
        sp_cols = [col for col in df_row.columns if sp in col]            
        df_y['theta_k '+sp] = np.vstack([fill,df_row[sp_cols].T.to_numpy()])
        
    # Nafion phase species' densities       
    for sp in ca.naf_b[0].species_names:
        sp_cols= [col for col in df_row.columns if sp in col]
        
        for j in range(cl.d['Nr']):
            r_cols = [col for col in sp_cols if 'r'+str(j) in col]            
            df_y['rho_k '+'r'+str(j)+' '+sp] = np.vstack([fill,df_row[r_cols].T.to_numpy()]) 
        
    # Faradaic current fraction
    if i_find != 0:
        sv = df_i.loc[df_ind].to_numpy()[1:]
    
        mvr = 0
        i_far_frac = np.zeros([cl.d['Ny'],1])
        
        for i in range(cl.d['Ny']):
            ca.inner_rxn_state(cl,sv,i)
            
            i_far_frac[i] = -ca.pt_s[i].get_net_production_rates(ca.carb)*ct.faraday\
                          *cl.d['SApv_pt'][i] / (df_row.iloc[0,0] *cl.d['1/dy'] *100**2)
            
            mvr = mvr +cl.d['nxt_y']
            
        df_y['i_far_frac [-]'] = np.vstack([fill,i_far_frac]) 
    
    """ Data frame to store steady-state radial gradients: """
    df_r = pd.DataFrame()
    
    for i in range(cl.d['Ny']):
        df_r['Radius y'+str(i)+' [nm]'] = 1/cl.d['1/r_j'][i,:] *1e9
    
    for sp in ca.naf_b[0].species_names:
        sp_cols = [col for col in df_row.columns if sp in col]
        
        for j in range(cl.d['Ny']):
            y_cols = [col for col in sp_cols if 'y'+str(j) in col]
            df_r['rho_k '+'y'+str(j)+' '+sp] = df_row[y_cols].T.to_numpy()  
    
    return df_f,df_p,df_y,df_r
