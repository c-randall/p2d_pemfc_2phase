""" Import needed modules """
"-----------------------------------------------------------------------------"
import numpy as np
import cantera as ct
from math import exp, pi
from scipy.optimize import fsolve
from Shared_Funcs.pem_pak import *
from scipy.sparse import coo_matrix
from Shared_Funcs.pemfc_property_funcs import *  
    
""" Pre-load Phases and Set States """
"-----------------------------------------------------------------------------"
an = electrode(name='anode')
ca = electrode(name='cathode')

# yaml file for common species/phases
yamlfile = 'Inputs/pemfc_cs_5nm.yaml'

# Anode Phases:
an.add_phase('carb',ct.Solution(yamlfile,'metal'))
an.carb.TP = T_an,P_an

an.add_phase('gas',ct.Solution(yamlfile,'anode_gas'))
an.gas.TP = T_an,P_an

an.add_phase('naf_b',ct.Solution(yamlfile,'naf_bulk_an'))
an.naf_b.TP = T_an,P_an

sp = [an.naf_b,an.gas]
an.add_phase('naf_s',ct.Interface(yamlfile,'naf_surf_an',sp))
an.naf_s.TP = T_an,P_an

sp = [an.carb,an.naf_b,an.gas]
an.add_phase('pt_s',ct.Interface(yamlfile,'Pt_surf_an',sp))
an.pt_s.TP = T_an,P_an

# Cathode Phases:
ca.add_phase('carb',ct.Solution(yamlfile,'metal'))
ca.carb.TP = T_ca,P_ca
    
ca.add_phase('gas',ct.Solution(yamlfile,'cathode_gas'))
ca.gas.TP = T_ca,P_ca

ca.add_phase('h2o_b',ct.Solution(yamlfile,'h2o_bulk_ca'))
ca.h2o_b.TP = T_ca,P_ca

sp = [ca.h2o_b,ca.gas]
ca.add_phase('h2o_gs',ct.Interface(yamlfile,'h2o_gas_surf_ca',sp))
ca.h2o_gs.TP = T_ca,P_ca

t_yamls = np.array([5,7,12,18,20,42])
ca.var_phase(name='naf_b')
ca.var_phase(name='naf_gs')
ca.var_phase(name='naf_ws')
ca.var_phase(name='pt_s')

for i in range(Ny_cl):
    t_i = np.argmin(abs(t_yamls - t_naf[i]*1e9))
    y_i = 'Inputs/pemfc_cs_'+str(t_yamls[t_i])+'nm.yaml'
    
    ca.naf_b[i] = ct.Solution(y_i,'naf_bulk_ca')
    ca.naf_b[i].TP = T_ca,P_ca
    
    # H2O_n = ca.naf_b[i].species('H2O(Naf)')
    # fuel_coeffs = H2O_n.thermo.coeffs
    # fuel_coeffs[[6,13]] = thermo_variable
    # H2O_n.thermo = ct.NasaPoly2(273.15,600.0,ct.one_atm,fuel_coeffs)
    # ca.naf_b[i].modify_species(ca.naf_b[i].species_index('H2O(Naf)'),H2O_n)
    
    sp = [ca.naf_b[i],ca.gas]
    ca.naf_gs[i] = ct.Interface(y_i,'naf_gas_surf_ca',sp)
    ca.naf_gs[i].TP = T_ca,P_ca
    
    sp = [ca.naf_b[i],ca.h2o_b]
    ca.naf_ws[i] = ct.Interface(y_i,'naf_h2o_surf_ca',sp)
    ca.naf_ws[i].TP = T_ca,P_ca

    sp = [ca.carb,ca.naf_b[i],ca.gas]
    ca.pt_s[i] = ct.Interface(y_i,'Pt_surf_ca_2s',sp)
    ca.pt_s[i].TP = T_ca,P_ca  

# Change basis to ensure TDY density takes in mass units:
basis = 'mass'

# Set each phase to mass basis:
ca.carb.basis = basis
ca.gas.basis  = basis

for i in range(Ny_cl):
    ca.naf_b[i].basis = basis
    ca.naf_gs[i].basis = basis
    ca.naf_ws[i].basis = basis
    ca.pt_s[i].basis  = basis

""" Pointers for Generality """
"-----------------------------------------------------------------------------"
# Pointers in Solution Vector:
gas_k, naf_k, pt_k = ca.gas.n_species, ca.naf_b[0].n_species, ca.pt_s[0].n_species

" Gas diffusion layer (GDL) variables "
ca.set_ptr('eps_w_gdl',0)
ca.set_ptr('rho_gdl_k',np.arange(1,gas_k+1))

" Catalyst layer (CL) variables "
ca.set_ptr('phi_dl',(ca.ptr['rho_gdl_k'][-1]+1)*Ny_gdl)
ca.set_ptr('eps_w_cl',(ca.ptr['rho_gdl_k'][-1]+1)*Ny_gdl+1)
ca.set_ptr('rho_gas_k',np.arange(ca.ptr['eps_w_cl']+1,
                                 ca.ptr['eps_w_cl']+1+gas_k))
ca.set_ptr('theta_pt_k',np.arange(ca.ptr['rho_gas_k'][-1]+1,
                                  ca.ptr['rho_gas_k'][-1]+1+pt_k))
ca.set_ptr('rho_naf_k',np.arange(ca.ptr['theta_pt_k'][-1]+1,
                                 ca.ptr['theta_pt_k'][-1]+1+naf_k))   

""" Pre-Processing and Initialization """
"-----------------------------------------------------------------------------"
" Determine SV length/initialize "
# Initialize all nodes according to cti and user inputs:
ca.gas.TPX = T_BC, P_BC, X_BC
SV0_gas_k  = ca.gas.density_mass*ca.gas.Y
SV0_naf_k  = ca.naf_b[0].density_mass*ca.naf_b[0].Y
SV0_pt_k   = ca.pt_s[0].coverages

SV0_gdl = np.tile(np.hstack((1e-3,SV0_gas_k)),Ny_gdl)

SV0_r  = np.tile(SV0_naf_k,Nr_cl) 
SV0_cl = np.tile(np.hstack((phi_ca_init,1e-3,SV0_gas_k,SV0_pt_k,SV0_r)),Ny_cl)

SV0 = np.zeros(len(SV0_gdl) + len(SV0_cl))
SV0 = np.hstack([SV0_gdl,SV0_cl])

L_gdl = int(len(SV0_gdl))
L_cl  = int(len(SV0_cl))
L_sv  = L_gdl + L_cl

" Geometric parameters "
# Performs calculations using Pt-loading information and assumes 1/2 spheres of
# Pt on the carbon surface. 

# After using the max(t_naf) to calulate SApv_dl, which is not a function of
# depth, the varialb t_naf values are used to calculate the depth-dependent 
# geometric parameters.

""" Core-shell geometric options """
"-----------------------------------------------------------------------------"  
# if not t_naf_runner:
class weighted_pt():
    def __init__(self):
        self.d = {}
        self.d['wt'] = np.ones(Ny_cl) *1/Ny_cl
        self.d['x'] = np.linspace(1,Ny_cl,Ny_cl)
        
    def linear(self):
        m = 2 / (Ny_cl**2 + Ny_cl)
        self.d['wt'] = m*self.d['x']
        
    def exponential(self):
        
        def func_b(b,Ny):
            ans = 0
            for i in range(Ny):
                ans = ans + b**(Ny-i)
            return ans - 1
        
        b = fsolve(lambda b: func_b(b,Ny_cl), 0.5)
        self.d['wt'] = b**(Ny_cl -self.d['x'] +1)        
    
a_n = weighted_pt()

if pt_dist == 'lin': a_n.linear()
if pt_dist == 'exp': a_n.exponential()

geom_out = rxn_areas_cs(w_Pt,t_ref,eps_g_cl,t_naf,r_c,r_Pt,rho_Pt,theta,Ny_cl,a_n)

SApv_dl  = geom_out['SApv_dl']
SApv_pt  = geom_out['SApv_pt']    
SA_dl_part = geom_out['SA_dl_part']
N_pt_part = geom_out['N_pt_part']

p_Pt = SApv_pt / SApv_dl *100
p_eff_SAnaf = geom_out['p_eff_SAnaf']
    
eps_g_cl, eps_n_cl, SApv_n = np.zeros(Ny_cl), np.zeros(Ny_cl), np.zeros(Ny_cl)
        
for i in range(Ny_cl):
    V_s_part = 4/3*np.pi*(r_c+t_naf[i])**3 
    eps_g_cl[i] = 1 - (SApv_dl[i] *V_s_part/ SA_dl_part[i])
    eps_n_cl[i] = ((r_c+t_naf[i])**3 -r_c**3) *(1 -eps_g_cl[i]) / (r_c+t_naf[i])**3
    SApv_n[i] = 3*(1 -eps_g_cl[i]) / (r_c+t_naf[i])
    
# Print geometry information to console:        
print('\nt_naf:',np.round(t_naf*1e9,3),'\neps_g_cl:',np.round(eps_g_cl,3))
print('eps_n_cl:',np.round(eps_n_cl,3),'\nSApv_n:',np.round(SApv_n,3))
print('SApv_dl:',np.round(SApv_dl,3),'\nSApv_Pt:',np.round(SApv_pt,3))
print('%Pt:',np.round(p_Pt,3),'\n')

# Tortuosity calculation via Bruggeman correlation [-]:
tau_g_gdl = eps_g_gdl**(-0.5)
tau_g_cl = eps_g_cl**(-0.5)
tau_n_cl = eps_n_cl**(-0.5)               

# Radius vectors for diffusion calculations [m]
dr = t_naf/Nr_cl
r_j = np.zeros([Ny_cl,Nr_cl])
t_shl = np.zeros([Ny_cl,Nr_cl])

for i in range(Ny_cl):
    r_j[i,:] = np.linspace((r_c+t_naf[i])-dr[i]/2, r_c+dr[i]/2, Nr_cl)
    t_shl[i,:] = np.tile(t_naf[i]/Nr_cl, Nr_cl)

r_jph = np.zeros([Ny_cl,Nr_cl -1])  # jph -> "j plus half"
for i in range(Nr_cl -1):
    r_jph[:,i] = np.mean(r_j[:,i:i+2])
    
# Vol fracs of Nafion shells for production terms [-]:
V_shl_tot = (r_c +np.reshape(t_naf,(Ny_cl,1)))**3 - r_c**3
Vf_shl = ((r_j +t_shl/2)**3 - (r_j -t_shl/2)**3) / V_shl_tot

" Calculate the anode equilibrium potential for polarization curves "
dgibbs_an  = an.pt_s.delta_gibbs
dphi_eq_an = -dgibbs_an / (n_elec_an*ct.faraday)

" Load the GDL parameters into a dictionary "
gdl = domain(name='gdl')
gdl.d['wt1'] = 0.5
gdl.d['wt2'] = 0.5
gdl.d['y'] = t_gdl
gdl.d['Ny'] = Ny_gdl
gdl.d['Len'] = L_gdl
gdl.d['nxt_y'] = int(L_gdl / Ny_gdl)
gdl.d['TPX_BC'] = T_BC, P_BC, X_BC
gdl.d['T'] = T_ca
gdl.d['eps_g'] = eps_g_gdl*np.ones(Ny_gdl)
gdl.d['eps_go'] = eps_g_gdl*np.ones(Ny_gdl)
gdl.d['dy'] = t_gdl / Ny_gdl
gdl.d['1/dy'] = Ny_gdl / t_gdl
gdl.d['mu_w'] = mu_w
gdl.d['ca_w'] = ca_gdl
gdl.d['st_w'] = st_h2o
gdl.d['r_p'] = r_p_gdl
gdl.d['K'] = 6e-12 / 0.75 *gdl.d['eps_go']

" Load the CL parameters into a dictionary "
cl = domain(name='cl')
cl.d['wt1'] = 0.5
cl.d['wt2'] = 0.5
cl.d['dy'] = t_cl / Ny_cl
cl.d['1/dy'] = Ny_cl / t_cl
cl.d['eps_n'] = eps_n_cl
cl.d['1/eps_n'] = 1 / eps_n_cl
cl.d['eps_g'] = np.copy(eps_g_cl)
cl.d['eps_go'] = np.copy(eps_g_cl)
cl.d['eps/tau2_n'] = eps_n_cl / tau_n_cl**2 
cl.d['y'] = t_cl
cl.d['Ny'] = Ny_cl
cl.d['Nr'] = Nr_cl
cl.d['Len'] = L_cl
cl.d['nxt_y'] = int(L_cl / Ny_cl)          
cl.d['nxt_r'] = naf_k                   
cl.d['iH_n'] = ca.naf_b[0].species_index('H(Naf)')
cl.d['iO2_n'] = ca.naf_b[0].species_index('O2(Naf)')
cl.d['iH2O_n'] = ca.naf_b[0].species_index('H2O(Naf)')
cl.d['SApv_pt'] = SApv_pt
cl.d['SApv_n'] = SApv_n
cl.d['1/CA_dl'] = 1 / (C_dl*SApv_dl) 
cl.d['p_eff_SAnaf'] = p_eff_SAnaf
cl.d['1/r_j'] = 1 / r_j 
cl.d['r_jph'] = r_jph
cl.d['1/t_shl'] = 1 / t_shl
cl.d['1/dr'] = 1 / dr
cl.d['Vf_shl'] = Vf_shl
cl.d['1/Vf_shl'] = 1 / Vf_shl
cl.d['1/gamma'] = 1 / ca.pt_s[0].site_density
cl.d['w_Pt'] = w_Pt
cl.d['T'] = T_ca
cl.d['mu_w'] = mu_w
cl.d['ca_w'] = ca_cl
cl.d['st_w'] = st_h2o
cl.d['r_p'] = r_p_cl
cl.d['SA_tot_part'] = 4*np.pi*(r_c+t_naf)**2
cl.d['N_pt_part'] = N_pt_part
cl.d['SA_dl_part'] = SA_dl_part
cl.d['SApv_dl'] = SApv_dl
cl.d['t_naf'] = t_naf
cl.d['p_Pt'] = p_Pt
cl.d['sig_method'] = sig_meth
cl.d['D_O2_method'] = D_O2_meth
cl.d['r_c'] = r_c
cl.d['r_w_c'] = np.sqrt(SApv_n*SA_dl_part/SApv_dl/N_pt_part/pi) 
cl.d['K'] = 8e-16 / 0.4 *cl.d['eps_go']

cl.d['D_naf_k'] = np.zeros(naf_k)
cl.d['D_naf_k'][cl.d['iH_n']] = D_naf_k_func(cl.d)
cl.d['D_naf_k'][cl.d['iO2_n']] = D_naf_k_func(cl.d)
cl.d['D_naf_k'][cl.d['iH2O_n']] = D_naf_k_func(cl.d)
    
" Calculate BC parameters between GDL and CL "
gdl_cl = domain(name='gdl_cl')
gdl_cl.d['dy'] = 0.5*gdl.d['dy'] + 0.5*cl.d['dy']
gdl_cl.d['1/dy'] = 1 / gdl_cl.d['dy']
gdl_cl.d['wt1'] = 0.5*gdl.d['dy'] / gdl_cl.d['dy']
gdl_cl.d['wt2'] = 0.5*cl.d['dy'] / gdl_cl.d['dy']

" Porosity func to get vals based on V_w "
gdl.update(ca,SV0)
cl.update(ca,SV0)

" Load any PEM parameters into a dictionary "
pem = domain(name='pem')
pem.d['R_naf'] = R_naf

" Set Jacobian Sparsity for speed increases "
sparsity = np.zeros([SV0.size,SV0.size])

# GDL eps_w
r = ca.ptr['eps_w_gdl']
sparsity[r,r] = 1.

c = ca.ptr['eps_w_gdl'] +gdl.d['nxt_y']
sparsity[r,c] = 1.

c = ca.ptr['rho_gdl_k']
sparsity[r,c] = 1.

c = ca.ptr['rho_gdl_k'] +gdl.d['nxt_y']
sparsity[r,c] = 1.

for i in range(1,gdl.d['Ny']-1):
    r = ca.ptr['eps_w_gdl'] +i*gdl.d['nxt_y']
    sparsity[r,r] = 1.
    
    c = ca.ptr['eps_w_gdl'] +(i-1)*gdl.d['nxt_y']
    sparsity[r,c] = 1.
    
    c = ca.ptr['eps_w_gdl'] +(i+1)*gdl.d['nxt_y']
    sparsity[r,c] = 1.
    
    c = ca.ptr['rho_gdl_k'] +i*gdl.d['nxt_y']
    sparsity[r,c] = 1.
    
    c = ca.ptr['rho_gdl_k'] +(i-1)*gdl.d['nxt_y']
    sparsity[r,c] = 1.
    
    c = ca.ptr['rho_gdl_k'] +(i+1)*gdl.d['nxt_y']
    sparsity[r,c] = 1.
    
r = ca.ptr['eps_w_gdl'] +(gdl.d['Ny']-1)*gdl.d['nxt_y']
sparsity[r,r] = 1.

c = ca.ptr['eps_w_gdl'] +(gdl.d['Ny']-2)*gdl.d['nxt_y']
sparsity[r,c] = 1.

c = ca.ptr['eps_w_cl']
sparsity[r,c] = 1.

c = ca.ptr['rho_gdl_k'] +(gdl.d['Ny']-1)*gdl.d['nxt_y']
sparsity[r,c] = 1.

c = ca.ptr['rho_gdl_k'] +(gdl.d['Ny']-2)*gdl.d['nxt_y']
sparsity[r,c] = 1.

c = ca.ptr['rho_gas_k']
sparsity[r,c] = 1.
    
# GDL rho_g_k
for sp in ca.gas.species_names:
    
    sp_ind = ca.gas.species_index(sp)
    
    r = ca.ptr['rho_gdl_k'][sp_ind]
    sparsity[r,r] = 1.
    
    c = ca.ptr['rho_gdl_k']
    sparsity[r,c] = 1.
    
    c = ca.ptr['rho_gdl_k'] +gdl.d['nxt_y']
    sparsity[r,c] = 1.
    
    c = ca.ptr['eps_w_gdl']
    sparsity[r,c] = 1.

    for i in range(1,gdl.d['Ny']-1):
        r = ca.ptr['rho_gdl_k'][sp_ind] +i*gdl.d['nxt_y']
        sparsity[r,r] = 1.
        
        c = ca.ptr['rho_gdl_k'] +i*gdl.d['nxt_y']
        sparsity[r,c] = 1.
        
        c = ca.ptr['rho_gdl_k'] +(i-1)*gdl.d['nxt_y']
        sparsity[r,c] = 1.
        
        c = ca.ptr['rho_gdl_k'] +(i+1)*gdl.d['nxt_y']
        sparsity[r,c] = 1.
        
        c = ca.ptr['eps_w_gdl'] +i*gdl.d['nxt_y']
        sparsity[r,c] = 1.
        
    r = ca.ptr['rho_gdl_k'][sp_ind] +(gdl.d['Ny']-1)*gdl.d['nxt_y']
    sparsity[r,r] = 1.
    
    c = ca.ptr['rho_gdl_k'] +(gdl.d['Ny']-1)*gdl.d['nxt_y']
    sparsity[r,c] = 1.
    
    c = ca.ptr['rho_gdl_k'] +(gdl.d['Ny']-2)*gdl.d['nxt_y']
    sparsity[r,c] = 1.
    
    c = ca.ptr['rho_gas_k']
    sparsity[r,c] = 1.
    
    c = ca.ptr['eps_w_gdl'] +(gdl.d['Ny']-1)*gdl.d['nxt_y']
    sparsity[r,c] = 1.    
        
# CL phi_dl
r = ca.ptr['phi_dl']
sparsity[r,r] = 1.

c = ca.ptr['phi_dl'] +cl.d['nxt_y']
sparsity[r,c] = 1.

ih2o_n = ca.naf_b[0].species_index('H2O(Naf)')
c = [ca.ptr['rho_naf_k'][ih2o_n] +j*cl.d['nxt_r'] for j in range(cl.d['Nr'])]
sparsity[r,c] = 1.

ih2o_n = ca.naf_b[1].species_index('H2O(Naf)')
c = [ca.ptr['rho_naf_k'][ih2o_n] +cl.d['nxt_y'] +j*cl.d['nxt_r'] for j in range(cl.d['Nr'])]
sparsity[r,c] = 1.

c = ca.ptr['theta_pt_k']
sparsity[r,c] = 1.

c = ca.ptr['rho_naf_k'] +(cl.d['Nr']-1)*cl.d['nxt_r']
sparsity[r,c] = 1.

for i in range(1,cl.d['Ny']-1):
    r = ca.ptr['phi_dl'] +i*cl.d['nxt_y']
    sparsity[r,r] = 1.
    
    c = ca.ptr['phi_dl'] +(i-1)*cl.d['nxt_y']
    sparsity[r,c] = 1.
    
    c = ca.ptr['phi_dl'] +(i+1)*cl.d['nxt_y']
    sparsity[r,c] = 1.
    
    ih2o_n = ca.naf_b[i].species_index('H2O(Naf)')
    c = [ca.ptr['rho_naf_k'][ih2o_n] +i*cl.d['nxt_y'] +j*cl.d['nxt_r'] for j in range(cl.d['Nr'])]
    sparsity[r,c] = 1.
    
    ih2o_n = ca.naf_b[i+1].species_index('H2O(Naf)')
    c = [ca.ptr['rho_naf_k'][ih2o_n] +(i+1)*cl.d['nxt_y'] +j*cl.d['nxt_r'] for j in range(cl.d['Nr'])]
    sparsity[r,c] = 1.
    
    c = ca.ptr['theta_pt_k'] +i*cl.d['nxt_y']
    sparsity[r,c] = 1.
    
    c = ca.ptr['rho_naf_k'] +i*cl.d['nxt_y'] +(cl.d['Nr']-1)*cl.d['nxt_r']
    sparsity[r,c] = 1.
    
r = ca.ptr['phi_dl'] +(cl.d['Ny']-1)*cl.d['nxt_y']
sparsity[r,r] = 1.

c = ca.ptr['phi_dl'] +(cl.d['Ny']-2)*cl.d['nxt_y']
sparsity[r,c] = 1.

ih2o_n = ca.naf_b[cl.d['Ny']-1].species_index('H2O(Naf)')
c = [ca.ptr['rho_naf_k'][ih2o_n] +(cl.d['Ny']-1)*cl.d['nxt_y'] +j*cl.d['nxt_r'] for j in range(cl.d['Nr'])]
sparsity[r,c] = 1.

c = ca.ptr['theta_pt_k'] +(cl.d['Ny']-1)*cl.d['nxt_y']
sparsity[r,c] = 1.

c = ca.ptr['rho_naf_k'] +(cl.d['Ny']-1)*cl.d['nxt_y'] +(cl.d['Nr']-1)*cl.d['nxt_r']
sparsity[r,c] = 1.

# CL eps_w
r = ca.ptr['eps_w_cl']
sparsity[r,r] = 1.

c = ca.ptr['eps_w_gdl'] +(gdl.d['Ny']-1)*gdl.d['nxt_y']
sparsity[r,c] = 1.

c = ca.ptr['eps_w_cl'] +cl.d['nxt_y']
sparsity[r,c] = 1.

c = ca.ptr['rho_gas_k']
sparsity[r,c] = 1.

c = ca.ptr['rho_gdl_k'] +(gdl.d['Ny']-1)*gdl.d['nxt_y']
sparsity[r,c] = 1.

c = ca.ptr['rho_gas_k'] +cl.d['nxt_y']
sparsity[r,c] = 1.

c = ca.ptr['rho_naf_k']
sparsity[r,c] = 1.

for i in range(1,cl.d['Ny']-1):
    r = ca.ptr['eps_w_cl'] +i*cl.d['nxt_y']
    sparsity[r,r] = 1. 
    
    c = ca.ptr['eps_w_cl'] +(i-1)*cl.d['nxt_y']
    sparsity[r,c] = 1.
    
    c = ca.ptr['eps_w_cl'] +(i+1)*cl.d['nxt_y']
    sparsity[r,c] = 1.
    
    c = ca.ptr['rho_gas_k'] +i*cl.d['nxt_y']
    sparsity[r,c] = 1.
    
    c = ca.ptr['rho_gas_k'] +(i-1)*cl.d['nxt_y']
    sparsity[r,c] = 1.
    
    c = ca.ptr['rho_gas_k'] +(i+1)*cl.d['nxt_y']
    sparsity[r,c] = 1.
    
    c = ca.ptr['rho_naf_k'] +i*cl.d['nxt_y']
    sparsity[r,c] = 1.
    
r = ca.ptr['eps_w_cl'] +(cl.d['Ny']-1)*cl.d['nxt_y']
sparsity[r,r] = 1.

c = ca.ptr['eps_w_cl'] +(cl.d['Ny']-2)*cl.d['nxt_y']
sparsity[r,c] = 1.

c = ca.ptr['rho_gas_k'] +(cl.d['Ny']-1)*cl.d['nxt_y']
sparsity[r,c] = 1.

c = ca.ptr['rho_gas_k'] +(cl.d['Ny']-2)*cl.d['nxt_y']
sparsity[r,c] = 1.

c = ca.ptr['rho_naf_k'] +(cl.d['Ny']-1)*cl.d['nxt_y']
sparsity[r,c] = 1.

# CL rho_g_k
for sp in ca.gas.species_names:
    
    sp_ind = ca.gas.species_index(sp)
    
    r = ca.ptr['rho_gas_k'][sp_ind]
    sparsity[r,r] = 1.
    
    c = ca.ptr['rho_gas_k']
    sparsity[r,c] = 1.
    
    c = ca.ptr['rho_gdl_k'] +(gdl.d['Ny']-1)*gdl.d['nxt_y']
    sparsity[r,c] = 1.
    
    c = ca.ptr['rho_gas_k'] +cl.d['nxt_y']
    sparsity[r,c] = 1.
    
    c = ca.ptr['eps_w_cl']
    sparsity[r,c] = 1.
    
    c = ca.ptr['rho_naf_k']
    sparsity[r,c] = 1.
    
    for i in range(1,cl.d['Ny']-1):
        r = ca.ptr['rho_gas_k'][sp_ind] +i*cl.d['nxt_y']
        sparsity[r,r] = 1.
        
        c = ca.ptr['rho_gas_k'] +i*cl.d['nxt_y']
        sparsity[r,c] = 1.

        c = ca.ptr['rho_gas_k'] +(i-1)*cl.d['nxt_y']
        sparsity[r,c] = 1.
        
        c = ca.ptr['rho_gas_k'] +(i+1)*cl.d['nxt_y']
        sparsity[r,c] = 1.
        
        c = ca.ptr['eps_w_cl'] +i*cl.d['nxt_y']
        sparsity[r,c] = 1.
        
        c = ca.ptr['rho_naf_k'] +i*cl.d['nxt_y']
        sparsity[r,c] = 1.
        
    r = ca.ptr['rho_gas_k'][sp_ind] +(cl.d['Ny']-1)*cl.d['nxt_y']
    sparsity[r,r] = 1.
    
    c = ca.ptr['rho_gas_k'] +(cl.d['Ny']-1)*cl.d['nxt_y']
    sparsity[r,c] = 1.
    
    c = ca.ptr['rho_gas_k'] +(cl.d['Ny']-2)*cl.d['nxt_y']
    sparsity[r,c] = 1.
    
    c = ca.ptr['eps_w_cl'] +(cl.d['Ny']-1)*cl.d['nxt_y']
    sparsity[r,c] = 1.
    
    c = ca.ptr['rho_naf_k'] +(cl.d['Ny']-1)*cl.d['nxt_y']
    sparsity[r,c] = 1.
    
# CL theta_pt_k
for sp in ca.pt_s[0].species_names:
    
    for i in range(cl.d['Ny']):
        
        sp_ind = ca.pt_s[i].species_index(sp)
        
        r = ca.ptr['theta_pt_k'][sp_ind] +i*cl.d['nxt_y']
        sparsity[r,r] = 1.
        
        c = ca.ptr['phi_dl'] +i*cl.d['nxt_y']
        sparsity[r,c] = 1.
        
        c = ca.ptr['theta_pt_k'] +i*cl.d['nxt_y']
        sparsity[r,c] = 1.
        
        c = ca.ptr['rho_naf_k'] +i*cl.d['nxt_y'] +(cl.d['Nr']-1)*cl.d['nxt_r']
        sparsity[r,c] = 1.
        
# CL rho_n_k
for sp in ca.naf_b[0].species_names:
    
    for i in range(cl.d['Ny']):
        
        sp_ind = ca.naf_b[i].species_index(sp)
    
        r = ca.ptr['rho_naf_k'][sp_ind] +i*cl.d['nxt_y']
        sparsity[r,r] = 1.
        
        c = ca.ptr['rho_naf_k'][sp_ind] +i*cl.d['nxt_y'] +cl.d['nxt_r']
        sparsity[r,c] = 1.
        
        c = ca.ptr['rho_naf_k'] +i*cl.d['nxt_y']
        sparsity[r,c] = 1.
        
        c = ca.ptr['rho_gas_k'] +i*cl.d['nxt_y']
        sparsity[r,c] = 1.
        
        c = ca.ptr['eps_w_cl'] +i*cl.d['nxt_y']
        sparsity[r,c] = 1.
        
        if sp == 'H(Naf)':
            sparsity[r,:] = 0.
            sparsity[r,r] = 1.
        
        for j in range(1,cl.d['Nr']-1):
            r = ca.ptr['rho_naf_k'][sp_ind] +i*cl.d['nxt_y'] +j*cl.d['nxt_r']
            sparsity[r,r] = 1.
            
            c = ca.ptr['rho_naf_k'][sp_ind] +i*cl.d['nxt_y'] +(j-1)*cl.d['nxt_r']
            sparsity[r,c] = 1.
            
            c = ca.ptr['rho_naf_k'][sp_ind] +i*cl.d['nxt_y'] +(j+1)*cl.d['nxt_r']
            sparsity[r,c] = 1.
            
            c = ca.ptr['eps_w_cl'] +i*cl.d['nxt_y']
            sparsity[r,c] = 1.
            
            if sp == 'H(Naf)':
                sparsity[r,:] = 0.
                sparsity[r,r] = 1.
            
        r = ca.ptr['rho_naf_k'][sp_ind] +i*cl.d['nxt_y'] +(cl.d['Nr']-1)*cl.d['nxt_r']
        sparsity[r,r] = 1.
        
        c = ca.ptr['rho_naf_k'][sp_ind] +i*cl.d['nxt_y'] +(cl.d['Nr']-2)*cl.d['nxt_r']
        sparsity[r,c] = 1.
        
        c = ca.ptr['rho_naf_k'] +i*cl.d['nxt_y'] +(cl.d['Nr']-1)*cl.d['nxt_r']
        sparsity[r,c] = 1.
        
        c = ca.ptr['theta_pt_k'] +i*cl.d['nxt_y']
        sparsity[r,c] = 1.
        
        c = ca.ptr['phi_dl'] +i*cl.d['nxt_y']
        sparsity[r,c] = 1.
        
        c = ca.ptr['eps_w_cl'] +i*cl.d['nxt_y']
        sparsity[r,c] = 1.
        
        if sp == 'H(Naf)':
            sparsity[r,:] = 0.
            sparsity[r,r] = 1.
        
# Sparse matrix
sparsity = coo_matrix(sparsity)