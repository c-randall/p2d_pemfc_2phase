""" Import needed modules """
"-----------------------------------------------------------------------------"
import numpy as np
import cantera as ct
from Shared_Funcs.pem_pak import *
from Shared_Funcs.pemfc_property_funcs import *  
    
""" Pre-load Phases and Set States """
"-----------------------------------------------------------------------------"
an = electrode(name='anode')
ca = electrode(name='cathode')

# Anode Phases:
an.add_phase('gas',ct.Solution(ctifile,'anode_gas'))
an.gas.TP = T_an,P_an

an.add_phase('naf_b',ct.Solution(ctifile,'naf_bulk_an'))
an.naf_b.TP = T_an,P_an

sp = [an.naf_b,an.gas]
an.add_phase('naf_s',ct.Interface(ctifile,'naf_surf_an',sp))
an.naf_s.TP = T_an,P_an

an.add_phase('carb',ct.Solution(ctifile,'metal'))
an.carb.TP = T_an,P_an

sp = [an.carb,an.naf_b,an.gas]
an.add_phase('pt_s',ct.Interface(ctifile,'Pt_surf_an',sp))
an.pt_s.TP = T_an,P_an

# Cathode Phases:
ca.add_phase('gas',ct.Solution(ctifile,'cathode_gas'))
ca.gas.TP = T_ca,P_ca

ca.add_phase('naf_b',ct.Solution(ctifile,'naf_bulk_ca'))
ca.naf_b.TP = T_ca,P_ca

sp = [ca.naf_b,ca.gas]
ca.add_phase('naf_s',ct.Interface(ctifile,'naf_surf_ca',sp))
ca.naf_s.TP = T_ca,P_ca

ca.add_phase('carb',ct.Solution(ctifile,'metal'))
ca.carb.TP = T_ca,P_ca

rm,sp = rxn_mech,[ca.carb,ca.naf_b,ca.gas]
ca.add_phase('pt_s',ct.Interface(ctifile,'Pt_surf_ca_'+rm,sp))
ca.pt_s.TP = T_ca,P_ca

# Change basis to ensure TDY density takes in mass units:
basis = 'mass'

# Set each phase to mass basis:
ca.gas.basis   = basis
ca.naf_s.basis = basis
ca.naf_b.basis = basis
ca.carb.basis  = basis
ca.pt_s.basis  = basis

""" Pointers for Generality """
"-----------------------------------------------------------------------------"
# Pointers in Solution Vector:
gas_k, naf_k, pt_k = ca.gas.n_species, ca.naf_b.n_species, ca.pt_s.n_species

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
ca.gas.TPY = T_ca_BC, P_ca_BC, Y_ca_BC
SV0_gas_k  = ca.gas.density_mass*ca.gas.Y
SV0_naf_k  = ca.naf_b.density_mass*ca.naf_b.Y
SV0_pt_k   = ca.pt_s.coverages

SV0_gdl = np.tile(np.hstack((0,SV0_gas_k)),Ny_gdl)

SV0_r  = np.tile(SV0_naf_k,Nr_cl) 
SV0_cl = np.tile(np.hstack((phi_ca_init,0,SV0_gas_k,SV0_pt_k,SV0_r)),Ny_cl)

SV0 = np.zeros(len(SV0_gdl) + len(SV0_cl))
SV0 = np.hstack([SV0_gdl,SV0_cl])

L_gdl = int(len(SV0_gdl))
L_cl  = int(len(SV0_cl))
L_sv  = L_gdl + L_cl

" Geometric parameters "
# Calculate all of the surface areas and volume fractions needed. Different
# methods include area_calcs = 0 which assumes that the Pt is flat circles and
# performs calculations with the specified % Pt. Using area_cals = 1 performs
# calculations using Pt-loading information and assumes 1/2 spheres of Pt on 
# the carbon surface. 

# Special cases for area_calcs = 2 and 3 were made for the shell thickness 
# runner of the core-shell model. They assume that as the Nafion shell becomes
# thinner, that the additional volume is given to the gas or carbon phase
# respectively.

""" Core-shell geometric options """
"-----------------------------------------------------------------------------"
if area_calcs == 0:
    # Determine method for Nafion SA based on theta:
    A_naf_reg  = 4*np.pi*(r_c+t_naf)**2
    A_naf_edit = 4*np.pi*r_c**2*(p_Pt/100)*(1+np.tan(np.deg2rad(theta)))**2
    
    # Naf/gas interface per total volume [m^2 Naf-gas int / m^3 tot]
    if A_naf_reg < A_naf_edit:
        SApv_naf = 3*(1 -eps_g_cl) / (r_c+t_naf)
        p_eff_SAnaf = 1.0
    else:
        SApv_naf = (1 -eps_g_cl) *A_naf_edit / (4/3 *np.pi *(r_c+t_naf)**3)
        p_eff_SAnaf = A_naf_edit / A_naf_reg
        
    # Double layer and Pt SA per total volume [m^2 cathode / m^3 tot]
    SApv_dl = 3*(1 -eps_g_cl) *r_c**2 / (r_c+t_naf)**3
    SApv_pt = 3*(1 -eps_g_cl) *r_c**2 *(p_Pt/100) / (r_c+t_naf)**3
    
    # Volume fraction of nafion in CL [-]
    eps_n_cl = ((r_c+t_naf)**3 -r_c**3) *(1 -eps_g_cl) / (r_c+t_naf)**3
        
elif area_calcs == 1:
    geom_out = rxn_areas_cs(w_Pt, t_cl, eps_g_cl, t_naf, r_c, r_Pt, rho_Pt, theta)
    
    SApv_naf = geom_out['SApv_naf']
    SApv_dl  = geom_out['SApv_dl']
    SApv_pt  = geom_out['SApv_pt']    
    eps_n_cl = geom_out['eps_n_cl']
    
    p_Pt = SApv_pt / SApv_dl *100
    p_eff_SAnaf = geom_out['p_eff_SAnaf']
        
elif area_calcs == 2: # gas phase changes
    eps_g_cl = 1 - (SApv_dl *(r_c+t_naf)**3 / (3 *r_c**2))
    eps_n_cl = ((r_c+t_naf)**3 -r_c**3) *(1 -eps_g_cl) / (r_c+t_naf)**3
    SApv_naf = 3*(1 -eps_g_cl) / (r_c+t_naf)
        
elif area_calcs == 3: # carbon phase changes
    r_c = 3*(1 -eps_g_cl) / SApv_naf - t_naf
    p_Pt = SApv_pt *(r_c+t_naf)**3 / (3*(1 -eps_g_cl) *r_c**2) *100
    SApv_dl = SApv_pt / p_Pt
    eps_n_cl = ((r_c+t_naf)**3 -r_c**3) *(1 -eps_g_cl) / (r_c+t_naf)**3

# Print geometry information to console:        
print('\nt_naf:',np.round(t_naf*1e9,3), '\neps_g_cl:',np.round(eps_g_cl,3))
print('eps_n_cl:',np.round(eps_n_cl,3),'\nA_int:',np.round(SApv_naf,3))
print('A_dl:',np.round(SApv_dl,3), '\nA_Pt:',np.round(SApv_pt,3))
print('%Pt:',np.round(p_Pt,3),'\n\n')

# Tortuosity calculation via Bruggeman correlation [-]:
tau_g_gdl = eps_g_gdl**(-0.5)
tau_g_cl = eps_g_cl**(-0.5)
if rxn_mech == '1s': tau_n_cl = 1.   # don't use Bruggeman for tau_naf
elif rxn_mech == '2s': tau_n_cl = 1. # don't use Bruggeman for tau_naf

# Radius vectors for diffusion calculations [m]
r_j = np.linspace((r_c+t_naf)-t_naf/(2*Nr_cl), r_c+t_naf/(2*Nr_cl), Nr_cl)
t_shl = np.tile(t_naf/Nr_cl, Nr_cl)
dr = abs(np.diff(r_j))

r_jph = np.zeros(Nr_cl -1) # jph -> "j plus half"
for i in range(Nr_cl -1):
    r_jph[i] = np.mean(r_j[i:i+2])
    
# Vol fracs of Nafion shells for weighted flux terms [-]:
Vf_shl = np.zeros(Nr_cl)
for i in range(Nr_cl):
    Vf_shl[i] = ((r_j[i] +t_shl[i]/2)**3 - (r_j[i] -t_shl[i]/2)**3) \
              / ((r_c+t_naf)**3 - r_c**3)

" Calculate the anode equilibrium potential for polarization curves "
dgibbs_an  = an.pt_s.delta_gibbs
dphi_eq_an = -dgibbs_an / (n_elec_an*ct.faraday)

" Initialize dictionaries to pass parameters to functions "
gdl, cl, gdl_cl, pem, p = {}, {}, {}, {}, {}

" Load the GDL parameters into a dictionary "
gdl['wt1'] = 0.5
gdl['wt2'] = 0.5
gdl['eps_g'] = eps_g_gdl
gdl['dy'] = t_gdl / Ny_gdl
gdl['1/dy'] = 1 / gdl['dy']
gdl['y'] = t_gdl
gdl['Ny'] = Ny_gdl
gdl['nxt_y'] = int(L_gdl / Ny_gdl)    # spacing between adjacent y nodes in GDL
gdl['TPY_BC'] = T_ca_BC, P_ca_BC, Y_ca_BC
gdl['T'] = T_ca
gdl['name'] = 'gdl'

gdl = porosity_func([0],gdl,0)

" Load the CL parameters into a dictionary "
cl['wt1'] = 0.5
cl['wt2'] = 0.5
cl['dy'] = t_cl / Ny_cl
cl['1/dy'] = 1 / cl['dy']
cl['1/eps_n'] = 1 / eps_n_cl
cl['eps_g'] = eps_g_cl
cl['eps/tau2_n'] = eps_n_cl / tau_n_cl**2 # based on full cell Nafion vol frac
cl['y'] = t_cl
cl['Ny'] = Ny_cl
cl['Nr'] = Nr_cl
cl['nxt_y'] = int(L_cl / Ny_cl)       # spacing between adjacent y nodes in CL
cl['nxt_r'] = naf_k                   # spacing between adjacent r nodes in CL
cl['i_e'] = ind_e
cl['iH_n'] = ca.naf_b.species_index('H(Naf)')
cl['iO2_n'] = ca.naf_b.species_index('O2(Naf)')
cl['iH2O_n'] = ca.naf_b.species_index('H2O(Naf)')
cl['SApv_pt'] = SApv_pt
cl['SApv_naf'] = SApv_naf
cl['1/CA_dl'] = 1 / (C_dl*SApv_dl) 
cl['p_eff_SAnaf'] = p_eff_SAnaf
cl['1/r_j'] = 1 / r_j 
cl['r_jph'] = r_jph
cl['1/t_shl'] = 1 / t_shl
cl['1/dr'] = 1 / dr
cl['Vf_shl'] = Vf_shl
cl['1/Vf_shl'] = 1 / Vf_shl
cl['1/gamma'] = 1 / ca.pt_s.site_density
cl['T'] = T_ca
cl['name'] = 'cl'

cl['D_naf_k'] = np.zeros(naf_k)
cl['D_naf_k'][cl['iH_n']] = D_O2_func(T_ca, t_naf, p_Pt, cl, D_O2_method)
cl['D_naf_k'][cl['iO2_n']] = D_O2_func(T_ca, t_naf, p_Pt, cl, D_O2_method)
cl['D_naf_k'][cl['iH2O_n']] = D_O2_func(T_ca, t_naf, p_Pt, cl, D_O2_method)
cl['sig_naf_io'] = sig_naf_io_func(T_ca, t_naf, RH, p_Pt, cl, sig_method)

cl = porosity_func([0],cl,0)
    
" Calculate BC parameters between GDL and CL "
gdl_cl['dy'] = 0.5*gdl['dy'] + 0.5*cl['dy']
gdl_cl['1/dy'] = 1 / gdl_cl['dy']
gdl_cl['wt1'] = 0.5*gdl['dy'] / gdl_cl['dy']
gdl_cl['wt2'] = 0.5*cl['dy'] / gdl_cl['dy']
gdl_cl['name'] = 'gdl_cl'

gdl_cl = porosity_func([0],gdl_cl,0,gdl=gdl,cl=cl)

" Load any PEM parameters into a dictionary "
pem['R_naf'] = R_naf

" Combine all dictionaries "
p['gdl'] = gdl
p['gdl_cl'] = gdl_cl
p['cl'] = cl
p['pem'] = pem