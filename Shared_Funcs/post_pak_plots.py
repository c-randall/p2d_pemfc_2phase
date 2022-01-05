"""
Functions to plot helpful information from pemfc model. Dataframes referrenced
in the plotting routine are established in post_pak_dfs, or separately loaded.
"""

import numpy as np
import cantera as ct
import matplotlib.pyplot as plt
from Shared_Funcs.read_and_write import *
from Shared_Funcs.pemfc_property_funcs import *
from Shared_Funcs.pemfc_transport_funcs import *

def fig_starter(fig_num):
    while plt.fignum_exists(fig_num):
        fig_num = fig_num +1
        
    return fig_num

def debug_plts(df_t,ca,gdl,cl,tog,save,save_folder):
    # Plot solution vector variables vs time to check for divergence issues
    # when trying to debug. Tog controls high and low detail.
    
    fig_num = fig_starter(0)
    
    """ GDL water volume fractions """
    eps_cols_gdl = []
    eps_cols_gdl.extend([col for col in df_t.columns if 'eps_w_gdl' in col])
    
    plt.figure(fig_num)
    plt.plot(df_t['Time [s]'],df_t[eps_cols_gdl])
    
    plt.legend(eps_cols_gdl,loc='best')
    plt.ylabel('GDL Water Volume Frac [-]')
    plt.xlabel('Time, t [s]')
    plt.tight_layout()
    
    if save:
        save_plot(os.getcwd(),save_folder,'GDL_eps_w_v_Time.png')
        
    fig_num = fig_num +1
    
    """ GDL gas densities """
    gdl_gas_cols = [col for col in df_t.columns if 'rho_gdl_k' in col]
    
    for i in range(gdl.d['Ny']):
        y_cols = [col for col in gdl_gas_cols if 'y'+str(i) in col]
        
        plt.figure(fig_num)
        plt.plot(df_t['Time [s]'],df_t[y_cols])
        
        plt.title('GDL y-node='+str(i))
        plt.legend(ca.gas.species_names,loc='best')
        plt.ylabel(r'GDL Gas $\rho_k$ [kg/m$^3$]')
        plt.xlabel('Time, t [s]')
        plt.tight_layout()    
        
        if save:
            fig_name = 'GDL_gas_densities_v_Time_y'+str(i)+'.png'
            save_plot(os.getcwd(),save_folder,fig_name)
        
        fig_num = fig_num +1
        
        if tog == 1:
            break
    
    """ Double layer potential """
    phi_cols = []
    phi_cols.extend([col for col in df_t.columns if 'phi_dl' in col])

    plt.figure(fig_num)
    plt.plot(df_t['Time [s]'],df_t[phi_cols])

    plt.legend(phi_cols,loc='best')
    plt.ylabel('Cathode DL Potential [V]')
    plt.xlabel('Time, t [s]')
    plt.tight_layout()
    
    if save:
        save_plot(os.getcwd(),save_folder,'Double_Layer_v_Time.png')
        
    fig_num = fig_num +1   
    
    """ CL water volume fractions """
    eps_cols_cl = []
    eps_cols_cl = [col for col in df_t.columns if 'eps_w_cl' in col]
    
    plt.figure(fig_num)
    plt.plot(df_t['Time [s]'],df_t[eps_cols_cl])
    
    plt.legend(eps_cols_cl,loc='best')
    plt.ylabel('CL Water Volume Frac [-]')
    plt.xlabel('Time, t [s]')
    plt.tight_layout()
    
    if save:
        save_plot(os.getcwd(),save_folder,'CL_eps_w_v_Time.png')
        
    fig_num = fig_num +1
        
    """ CL gas densities """
    cl_gas_cols = [col for col in df_t.columns if 'rho_gas_k' in col]
    
    for i in range(cl.d['Ny']):
        y_cols = [col for col in cl_gas_cols if 'y'+str(i) in col]
        
        plt.figure(fig_num)
        plt.plot(df_t['Time [s]'],df_t[y_cols])
        
        plt.title('CL y-node='+str(i))
        plt.legend(ca.gas.species_names,loc='best')
        plt.ylabel(r'CL Gas $\rho_k$ [kg/m$^3$]')
        plt.xlabel('Time, t [s]')
        plt.tight_layout() 
        
        if save:
            fig_name = 'CL_gas_densities_v_Time_y'+str(i)+'.png'
            save_plot(os.getcwd(),save_folder,fig_name)
        
        fig_num = fig_num +1
        
        if tog == 1:
            break
        
    """ CL Pt surface sites """
    cl_pt_cols = [col for col in df_t.columns if 'theta_pt_k' in col]
    
    for i in range(cl.d['Ny']):
        y_cols = [col for col in cl_pt_cols if 'y'+str(i) in col]
        
        plt.figure(fig_num)
        plt.plot(df_t['Time [s]'],df_t[y_cols])
        
        plt.title('CL y-node='+str(i))
        plt.legend(ca.pt_s[0].species_names)
        plt.ylabel('Surface Coverage [-]')
        plt.xlabel('Time, t [s]')
        plt.tight_layout()
        
        if save:
            fig_name = 'CL_pt_coverages_v_Time_y'+str(i)+'.png'
            save_plot(os.getcwd(),save_folder,fig_name)
        
        fig_num = fig_num +1
        
        if tog == 1:
            break      
        
    """ CL Nafion species' densities """
    cl_naf_cols = [col for col in df_t.columns if 'rho_naf_k' in col]

    for i in range(cl.d['Ny']):
        y_cols = [col for col in cl_naf_cols if 'y'+str(i) in col]
        
        for sp in ca.naf_b[0].species_names:            
            if sp == 'H(Naf)' and tog == 1:
                None
            else:
                sp_cols = [col for col in y_cols if sp in col]
                plt.figure(fig_num)
                
                for j in range(cl.d['Nr']):
                    plt_col = []
                    plt_col.extend([col for col in sp_cols if 'r'+str(j) in col])
                    plt.plot(df_t['Time [s]'],df_t[plt_col])
                    
                plt.title('CL y-node='+str(i))
                plt.legend(['r'+str(n) for n in range(cl.d['Nr'])],loc='best')
                plt.ylabel(sp+r' $\rho_k$ [kg/m$^3$]')
                plt.xlabel('Time, t [s]')
                plt.tight_layout()
        
                if save:
                    fig_name = 'CL_naf_densities_v_Time_y'+str(i)+'.png'
                    save_plot(os.getcwd(),save_folder,fig_name)
                
                fig_num = fig_num +1
            
        if tog == 1:
            break
        
    return None
        
def grad_plts(df_y,ca,gdl,cl,tog,save,save_folder):
    # Plot solution variable depth gradients to identify trends and limiting
    # transport phenomena. Tog controls high and low detail.
    
    fig_num = fig_starter(0)
    
    """ Water volume fraction """
    plt.figure(fig_num)
    eps_w = df_y['eps_w [-]'].to_numpy()
    s_w_gdl = eps_w[:gdl.d['Ny']] / gdl.d['eps_go']
    s_w_cl = eps_w[gdl.d['Ny']:] / cl.d['eps_go']
    plt.plot(df_y['Depth [um]'],np.hstack([s_w_gdl,s_w_cl]),'-o')
    
    plt.xlabel(r'Cathode Depth [$\mu$m]')
    plt.ylabel(r'Water saturation [-]')
    plt.tight_layout()
    
    if save:
        save_plot(os.getcwd(),save_folder,'Water_vol_frac_v_Depth.png')
    
    fig_num = fig_num +1
    
    """ Gas phase species' densities """
    for sp in ca.gas.species_names:
        plt.figure(fig_num)
        gas_k_col = [col for col in df_y.columns if sp+'(gas)' in col]
        plt.plot(df_y['Depth [um]'],df_y[gas_k_col],'-o')
        
        if tog == 2:
            plt.xlabel(r'Cathode Depth [$\mu$m]')
            plt.ylabel(sp+r'(gas) density [kg/m$^3$]')
            plt.tight_layout()
            
            if save:
                save_plot(os.getcwd(),save_folder,'Gas_Phase_'+sp+'_Density_v_Depth.png')
                
            fig_num = fig_num +1
                
    if tog == 1:
        plt.legend(ca.gas.species_names,loc='best')
        plt.xlabel(r'Cathode Depth [$\mu$m]')
        plt.ylabel(r'Gas phase $\rho_k$ [kg/m$^3$]')
        plt.tight_layout()

        if save:
            save_plot(os.getcwd(),save_folder,'Gas_Phase_Densities_v_Depth.png')
        
        fig_num = fig_num +1
    
    """ Double layer potential """
    plt.figure(fig_num)
    x_vals = (df_y['Depth [um]'][gdl.d['Ny']:] - gdl.d['y']*1e6)
    y_vals = -1*df_y['phi_dl'][gdl.d['Ny']:]
    plt.plot(x_vals,y_vals,'-o')
    
    plt.xlabel(r'Cathode CL Depth [$\mu$m]')
    plt.ylabel(r'Electrolyte Potential [V]')
    plt.tight_layout()
    
    if save:
        save_plot(os.getcwd(),save_folder,'Nafion_Potential_v_CL_Depth.png')
        
    fig_num = fig_num +1
    
    """ Pt surface coverages """
    for sp in ca.pt_s[0].species_names:
        plt.figure(fig_num)      
        
        sp_cols = []
        sp_cols.extend([col for col in df_y.columns if sp in col])
        
        x_vals = (df_y['Depth [um]'][gdl.d['Ny']:] - gdl.d['y']*1e6)
        y_vals = df_y[sp_cols][gdl.d['Ny']:]
        plt.plot(x_vals,y_vals,'-o')
        
        if tog == 2:
            plt.xlabel(r'Cathode CL Depth [$/mu$m]')
            plt.ylabel(sp+r' coverage [-]')
            plt.tight_layout()
            
            if save:
                save_plot(os.getcwd(),save_folder,sp+'_Coverage_v_CL_Depth.png')
                
            fig_num = fig_num +1  
    
    if tog == 1:
        plt.legend(ca.pt_s[0].species_names,loc='best')
        plt.xlabel(r'Cathode CL Depth [$/mu$m]')
        plt.ylabel(r'Surface coverage [-]')
        plt.tight_layout()
        
        if save:
            save_plot(os.getcwd(),save_folder,'Pt_Coverage_v_CL_Depth.png')
            
        fig_num = fig_num +1  
    
    """ Nafion phase species' densities """
    for sp in ca.naf_b[0].species_names:               
        if sp == 'H(Naf)' and tog == 1:
            None
        else:
            sp_cols = []
            sp_cols.extend([col for col in df_y.columns if sp in col])
        
            plt.figure(fig_num)
            
            for i in range(cl.d['Nr']):                
                sp_col_r = [col for col in sp_cols if 'r'+str(i) in col]
                x_vals = (df_y['Depth [um]'][gdl.d['Ny']:] - gdl.d['y']*1e6)
                y_vals = df_y[sp_col_r][gdl.d['Ny']:]
                
                plt.plot(x_vals,y_vals,'-o')
            
            plt.xlabel(r'Cathode CL Depth [$/mu$m]')
            plt.ylabel(sp+r' density [kg/m$^3$]')
            plt.legend(['r'+str(i) for i in range(cl.d['Nr'])])
            plt.tight_layout()
            
            if save:
                fig_name = 'Nafion_Phase_'+sp+'_Density_v_CL_Depth.png'
                save_plot(os.getcwd(),save_folder,fig_name)
            
            fig_num = fig_num +1            
        
    """ Faradaic current fraction """
    if 'i_far_frac [-]' in df_y.columns:
        plt.figure(fig_num)
        x_vals = (df_y['Depth [um]'][gdl.d['Ny']:] - gdl.d['y']*1e6)
        y_vals = df_y['i_far_frac [-]'][gdl.d['Ny']:]
        plt.plot(x_vals,y_vals,'-o')
        
        plt.xlabel(r'Cathode CL Depth [$/mu$m]')
        plt.ylabel(r'i$_{Far}$ / i$_{ext}$ [-]')
        plt.tight_layout()
            
        if save:
            save_plot(os.getcwd(),save_folder,'i_far_frac_v_CL_Depth.png')
            
        fig_num = fig_num +1
    
    return None
    
def radial_plts(df_r,ca,gdl,cl,tog,save,save_folder):
    # Plot solution variable radial gradients to identify trends and limiting
    # transport phenomena. Tog controls high and low detail.
    
    fig_num = fig_starter(0)
    cnt = 0
    
    for sp in ca.naf_b[0].species_names:        
        if sp == 'H(Naf)' and tog == 1:
            None
        else:
            sp_cols = []
            sp_cols.extend([col for col in df_r.columns if sp in col])
            
            for i in range(cl.d['Ny']):
                plt.figure(fig_num +cnt)
                
                plt_col_x = [col for col in df_r.columns if 'Radius y'+str(i) in col]
                plt_col_y = [col for col in sp_cols if 'y'+str(i) in col]
                plt.plot(df_r[plt_col_x],df_r[plt_col_y],'-o')
                
                if tog == 1:
                    break
                
                if tog == 3:                
                    plt.legend(['y'+str(i)],loc='best')
                    plt.xlabel(r'Nafion Shell Radius [nm]')
                    plt.ylabel(sp+r' density [kg/m$^3$]')
                    plt.tight_layout()
                    
                    cnt = cnt +1
                    
            if tog != 3:
                plt.legend(['y'+str(j) for j in range(cl.d['Ny'])],loc='best')
                plt.xlabel(r'Nafion Shell Radius [nm]')
                plt.ylabel(sp+r' density [kg/m$^3$]')
                plt.tight_layout()
            
            if save:
                save_plot(os.getcwd(),save_folder,sp+'Density_v_Shell_Radius.png')
        
            fig_num = fig_num +1
            
    return None
    
def polar_plt(df_p,ca,gdl,cl,polar,data,save,save_folder):
    # Plot the polarization curve. Tog controls high and low detail. If data
    # requested, overlay Owejan et al. after checking that w_Pt matches.
    
    fig_num = fig_starter(0)
    
    fig = plt.figure(fig_num)
    ax1 = fig.add_axes([0.18, 0.2, 0.7, 0.7])
    ax1.plot(df_p['i_ext [A/cm2]'],df_p['Voltage [V]'],label='model')

    ax1.set_ylabel(r'Cell Voltage [V]')
    ax1.set_xlabel(r'Current Density [A/cm$^2$]')
    
    if all([data == 'air', cl.d['w_Pt'] == 0.2]):
        y = np.array([0.95,0.85,0.80,0.77,0.73,0.72,0.70,0.68,0.67,0.65,0.63]) 
        s = np.array([0.1,12,7,7,12,1,8,7,7,9,9]) *1e-3 
    elif all([data == 'air', cl.d['w_Pt'] == 0.1]):
        y = np.array([0.93,0.83,0.79,0.75,0.71,0.69,0.67,0.65,0.64,0.62,0.60]) 
        s = np.array([0.1,9,7,5,7,11,11,7,9,11,11]) *1e-3
    elif all([data == 'air', cl.d['w_Pt'] == 0.05]):
        y = np.array([0.92,0.81,0.76,0.72,0.67,0.65,0.63,0.60,0.59,0.56,0.54]) 
        s = np.array([0.1,8,6,6,7,7,5,5,6,7,7]) *1e-3
    elif all([data == 'air', cl.d['w_Pt'] == 0.025]):
        y = np.array([0.91,0.79,0.72,0.68,0.63,0.60,0.57,0.53,0.50,0.46,0.43]) 
        s = np.array([0.1,4,10,14,13,13,19,24,25,23,24]) *1e-3
    elif all([data == 'o2', cl.d['w_Pt'] == 0.2]):
        y = np.array([0.90,0.89,0.86,0.84,0.81,0.77,0.76,0.75,0.73,0.71,0.70])
        s = np.array([5,4,8,5,5,5,5,5,6,9,10]) *1e-3
    elif all([data == 'o2', cl.d['w_Pt'] == 0.1]):
        y = np.array([0.89,0.86,0.84,0.81,0.78,0.74,0.73,0.71,0.69,0.68,0.67])
        s = np.array([5,9,5,5,5,5,5,5,8,9,10]) *1e-3
    elif all([data == 'o2', cl.d['w_Pt'] == 0.05]):
        y = np.array([0.86,0.83,0.81,0.78,0.75,0.71,0.69,0.67,0.65,0.63,0.62])
        s = np.array([8,8,6,6,7,8,8,8,9,8,7]) *1e-3
    elif all([data == 'o2', cl.d['w_Pt'] == 0.025]):
        y = np.array([0.84,0.81,0.78,0.75,0.72,0.67,0.65,0.64,0.61,0.59,0.57])
        s = np.array([6,5,5,5,8,12,13,14,16,18,20]) *1e-3
    else:
        x = y = s = None
        
    if data == 'air':
        x = np.array([0.0,0.05,0.20,0.40,0.80,1.0,1.2,1.5,1.65,1.85,2.0])
    elif data == 'o2':
        x = np.array([0.03,0.05,0.10,0.2,0.4,0.8,1.0,1.2,1.5,1.75,2.0])
        
    if data:
        ax1.errorbar(x,y,yerr=s,fmt='.',color='C0',capsize=3,label='Owejan')
        ax1.set_ylim([0.35, 1.0])
        ax1.set_xlim([0, 2.1])
        ax1.legend(loc='lower center')
        
    if polar == 2:
        ax2 = ax1.twinx()
        ax2.plot(df_p['i_ext [A/cm2]'],df_p['Power [W/cm2]'],'--',color='C0')
        ax2.set_ylabel(r'Power Density [W/cm$^2$]')
    
    if save:
        save_plot(os.getcwd(),save_folder,'Polarization_Curve.png')
    
    fig_num = fig_num +1
    
    return None

def over_plt(df_p,ca,gdl,cl,save,save_folder):
    # Plot the overpotential curve. 
    
    fig_num = fig_starter(0)
    
    plt.figure(fig_num)
    plt.plot(df_p['i_ext [A/cm2]'], df_p['Eta [V]'])

    plt.ylabel(r'Overpotential [V]')
    plt.xlabel(r'Current Density [A/cm$^2$]')
    plt.tight_layout()
    
    if save:
        save_plot(os.getcwd(),save_folder,'Overpotential_Curve.png')
    
    fig_num = fig_num +1
    
    return None

def verification(df_i,ca,gdl,cl,gdl_cl,i_find):
    i_ind = np.argmin(abs(df_i[0] - i_find))
    sv = df_i.loc[i_ind][1:].to_numpy()
    
    cl.update(ca,sv)
    gdl.update(ca,sv)

    i_4F = df_i[0][i_ind]*100**2 / (4*ct.faraday)    

    i_Last_gdl = int(gdl.d['Len'] / gdl.d['Ny'] *(gdl.d['Ny'] -1))
    rho_gdl_k = sv[ca.ptr['rho_gdl_k'] +i_Last_gdl]
    TDY1 = gdl.d['T'], sum(rho_gdl_k), rho_gdl_k

    i_First_cl = 0
    rho_cl_k = sv[ca.ptr['rho_gas_k'] +i_First_cl]
    TDY2 = cl.d['T'], sum(rho_cl_k), rho_cl_k

    TDY_vec = [TDY1, TDY2]
    O2_BC_flux = fickian_adf(TDY_vec, ca, gdl_cl, gdl.d, cl.d, None) \
              / ca.gas.molecular_weights

    print('\ni_ext:', np.round(df_i[0][i_ind],3))
    print('O2_i_4F(x10^5):', np.round(i_4F*1e5,3))
    print('O2_BC_flux(x10^5):', np.round(O2_BC_flux[0]*1e5,3))
    print('ratio:', np.round(i_4F / O2_BC_flux[0],3),'\n')   
    
    return None