"""
This tool is to be used with pemfc models. It generalizes the use of building
phases, pointers, geometry, and transport.
"""

import numpy as np
from math import pi, exp, tanh
from Shared_Funcs.pemfc_property_funcs import *  

# generic class to define electrodes
class electrode():
    
    def __init__(self,**kwargs):
        self.name = kwargs.get('name','electrode')
        self.ptr  = {}
        
    def var_phase(self,**kwargs):
        setattr(self,kwargs.get('name','phase'),{})
        
    def add_phase(self,phase_name,ct_obj):
        setattr(self,phase_name,ct_obj)
    
    def set_ptr(self,ptr_name,ptr_val):
        self.ptr[ptr_name] = ptr_val
    
    def outter_rxn_state(self,cl,sv,i):
        rho_gas_k = sv[self.ptr['rho_gas_k'] +i*cl.d['nxt_y']]
        self.gas.TDY = cl.d['T'], sum(rho_gas_k), rho_gas_k
        
        rho_naf_k = sv[self.ptr['rho_naf_k'] +i*cl.d['nxt_y']]
        self.naf_b[i].Y = rho_naf_k
    
    def inner_rxn_state(self,cl,sv,i):
        phi_ref, phi_elyte = 0., -sv[self.ptr['phi_dl'] +i*cl.d['nxt_y']]
        
        self.carb.electric_potential = phi_ref
        self.pt_s[i].electric_potential = phi_ref

        self.naf_b[i].electric_potential = phi_elyte
        self.naf_gs[i].electric_potential = phi_elyte
        self.naf_ws[i].electric_potential = phi_elyte
        
        theta_pt_k = sv[self.ptr['theta_pt_k'] +i*cl.d['nxt_y']]
        self.pt_s[i].coverages = theta_pt_k
        
        inner_index = i*cl.d['nxt_y'] +(cl.d['Nr']-1)*cl.d['nxt_r']
        rho_naf_k = sv[self.ptr['rho_naf_k'] + inner_index]
        self.naf_b[i].Y = rho_naf_k        
        
# variable domain geometry/transport
class domain():
    
    def __init__(self,**kwargs):
        self.name    = kwargs.get('name',None)
        self.cutoff  = 1e-15
        self.dropoff = 1e-10
        self.d = {}
    
    def update(self,ca,sv):
        if self.name == 'gdl':
            ptr = ca.ptr['eps_w_gdl']
        elif self.name == 'cl':
            ptr = ca.ptr['eps_w_cl']
            
        eps_w = np.zeros(self.d['Ny'])
        for i in range(self.d['Ny']):
            eps_w[i] = max(sv[ptr +i*self.d['nxt_y']], self.cutoff)
            
        self.d['eps_w'] = eps_w
        self.d['eps_g'] = self.d['eps_go'] - self.d['eps_w']
        self.d['s_w'] = 1 - self.d['eps_g'] / self.d['eps_go']
                        
        self.geom()
        self.trans(ca,sv)
       
    def geom(self):
        self.d['K_w'] = self.d['K'] *self.d['s_w']**3
        self.d['K_g'] = self.d['K'] *(1 - self.d['s_w'])**3
            
        if self.name == 'cl':
            
            self.d['SApv_wg'] = self.d['SApv_n'] *self.d['s_w']
            self.d['SApv_nw'] = self.d['SApv_n'] *self.d['s_w']
            self.d['SApv_ng'] = self.d['SApv_n'] *self.d['p_eff_SAnaf'] *(1 - self.d['s_w']) 
            
            # num = self.d['eps_w'] / self.d['SApv_dl'] *self.d['SA_dl_part']
            # den1 = 2/3 *self.d['N_pt_part'] *pi
            # den2 = 4/3 *pi
            
            # r_w_nuc = (num/den1)**(1/3)
            # r_w_shl = (num/den2)**(1/3) - self.d['r_c'] - self.d['t_naf']
            
            # self.d['r_w'] = r_w_nuc
            
            # ind = self.d['r_w'] >= self.d['r_w_c']
            # self.d['r_w'][ind] = r_w_shl[ind]
            
            # nw_int_nuc = self.d['N_pt_part']*pi*self.d['r_w']**2/self.d['SA_dl_part']*self.d['SApv_dl']
            # nw_int_shl = self.d['SApv_n']
            
            # self.d['SApv_nw'] = nw_int_nuc
            # self.d['SApv_nw'][ind] = nw_int_shl[ind]
            
            # wg_int_nuc = self.d['N_pt_part']*2*pi*self.d['r_w']**2/self.d['SA_dl_part']*self.d['SApv_dl']
            # wg_int_shl = 4*pi*(self.d['r_c']+self.d['t_naf']+self.d['r_w'])**2/self.d['SA_dl_part']*self.d['SApv_dl']
            
            # self.d['SApv_wg'] = wg_int_nuc
            # self.d['SApv_wg'][ind] = wg_int_shl[ind]
            
            # self.d['SApv_ng'] = self.d['SApv_n'] - self.d['SApv_nw']
            
        self.d['1/eps_g'] = 1 / self.d['eps_g']
        self.d['eps/tau2'] = self.d['eps_g']**2
    
    def trans(self,ca,sv):
        if self.name == 'cl':
            self.d['sig_io'] = sig_io_func(self.d, ca, sv)
    