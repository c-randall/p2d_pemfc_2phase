import numpy as np

"""
Fickian Diffusion function (w/ advective term) in 1-D cartesian:
    This function takes in the TDY states at two adjacent points and uses
    average conditions at the boundary between them to calculate the mass flux
    by combining a gradient in mass fractions (Fickian diffusion) with a
    pressure driven gradient (advection term). Hence the name adf for advection
    diffusion flux.

    Positive flux implies net flux from the TDY1 node into the TDY2 node. Make
    sure to correctly emphasize this sign convention in runner codes."""

def fickian_adf(TDY_vec, ca, p, p1, p2, i):
    " Inputs to this function are as follows: "
    # TDY_vec: temperature, density, and mass fractions [K, kg/m^3, -]
    # ca: cathode class containing cantera obect as 'gas' attribute
    # p: dictionary with 'wt1'(gdl), 'wt2'(cl), '1/dy'
    # p1: dictionary with 'K_g', 'eps/tau2' for node 1
    # p2: dictionary with 'K_g', 'eps/tau2' for node 2
    # i: number for node 1 to pull correct property indices
    
    " Descriptions of the dictionary terms are as follows: "
    # 'wt1', 'wt2': fractional dy of gdl and cl for weighted averages [-]
    # 'K_g': permeability of the porous domain [m^2]
    # '1/dy': inverse of distance between nodes [1/m]
    # 'eps/tau2': porosity of the domain divided by its tortuosity squared [-]
    
    TDY1, TDY2 = TDY_vec[0], TDY_vec[1]
    
    # set state 1 properties:
    ca.gas.TDY = TDY1
    D_k1 = ca.gas.mix_diff_coeffs_mass
    rho1 = ca.gas.density_mass
    mu1 = ca.gas.viscosity
    P1 = ca.gas.P
    Y1 = ca.gas.Y
    rho_k1 = rho1*Y1
    
    # set state 2 properties:
    ca.gas.TDY = TDY2
    D_k2 = ca.gas.mix_diff_coeffs_mass
    rho2 = ca.gas.density_mass
    mu2 = ca.gas.viscosity
    P2 = ca.gas.P
    Y2 = ca.gas.Y
    rho_k2 = rho2*Y2
    
    # calculate average boundary properties:
    " when not at a boundary, wt terms shoud each be 0.5 "
    
    if i == 'channel':
        half = 2
        K_g = p1['K_g'][0]
        epstau = p1['eps/tau2'][0]
    elif p.name == 'gdl_cl':
        half = 1
        K_g = p.d['wt1']*p1['K_g'][-1] + p.d['wt2']*p2['K_g'][0]
        epstau = p.d['wt1']*p1['eps/tau2'][-1] + p.d['wt2']*p2['eps/tau2'][0] 
    else:
        half = 1
        K_g = p.d['wt1']*p1['K_g'][i] + p.d['wt2']*p2['K_g'][i+1]
        epstau = p.d['wt1']*p1['eps/tau2'][i] + p.d['wt2']*p2['eps/tau2'][i+1]   
    
    D_k = p.d['wt1']*D_k1 + p.d['wt2']*D_k2
    rho = p.d['wt1']*rho1 + p.d['wt2']*rho2
    mu = p.d['wt1']*mu1 + p.d['wt2']*mu2
    rho_k = p.d['wt1']*rho_k1 + p.d['wt2']*rho_k2    
    
    # convective and diffusive driving terms:
    J_conv = -rho_k*K_g*(P2 - P1)*half*p.d['1/dy'] / mu
    J_diff = -epstau*D_k*rho*(Y2 - Y1)*half*p.d['1/dy']
    
    # net mass flux of each species:
    mass_flux = J_conv + J_diff
    
    # Output returns mass flux vector [kg/m^2-s]
    return mass_flux

"""
Fickian Diffusion function in 1-D radial:
    This function uses a finite difference method to approximate the radial 
    mass flux between two adjacent nodes. 
    
    When using this code:
    Positive flux implies net flux from the outer (rho_k1) node into the inner
    (rho_k2) node. Make sure to correctly emphasize this sign convention in
    runner codes."""
    
def radial_fdiff(rho_k1, rho_k2, p, ind):
    " Inputs to this function are as follows: "
    # rho_k1, rho_k1: vectors with mass density of each species [kg/m^3]
    # p: parameters with terms as defined below
    # ind: index [i,j] to pull correct geometric values from p
    
    " Descriptions of the dictionary terms are as follows: "
    # 'D_naf_k': effective duffusion coefficients [m^2/s]
    # 'r_jph': radius at bounday between nodes (j+1/2) [m]
    # '1/dr': inverse of distance between nodes [1/m]
    # 'p_eff_SAnaf': effective Nafion SA based on Pt and max angle [-]
    
    drho_dt = p.d['D_naf_k'] *p.d['p_eff_SAnaf'][ind[0]] *(1-p.d['s_w'][ind[0]])\
            *(p.d['r_jph'][ind[0],ind[1]]**2 *(rho_k1 -rho_k2) *p.d['1/dr'][ind[0]])
                	
    # Output returns change in mass vector [kg/s]
    return drho_dt

"""
Darcy's law for liquid water transport:
    This function uses a finite difference method to approximate the volume
    flux between two adjacent nodes in Cartesian coordinates.
    
    Positive flux implies net flux from the TDY1 node into the TDY2 node. Make
    sure to correctly emphasize this sign convention in runner codes."""
    
def J_hydrophilic(s_w):
    J = 1.417*(1-s_w) - 2.120*(1-s_w)**2 + 1.262*(1-s_w)**3
    return J

def J_hydrophobic(s_w):
    J = 1.417*s_w - 2.120*s_w**2 + 1.262*s_w**3
    return J

def P_cap(p):
    P_cap = p['st_w']*np.cos(np.deg2rad(p['ca_w']))*(p['eps_go']/p['K'])**0.5
    return P_cap
    
def darcys_law(TDY_vec, eps_vec, ca, p, p1, p2, i):
    " Inputs to this function are as follows: "
    # TDY_vec: temp, density, and mass fracs of gas phase [K, kg/m^3, -]
    # eps_vec: water volume fractions for node 1 and 2 [-]
    # ca: class holding cantera objects as 'gas' and 'h2o_b' attributes
    # p: dictionary with weighting terms for node 1 and node 2
    # p1: dictionary with 'K_w', 'P_c', 'mu_w', '1/dy' for node 1
    # p2: dictionary with 'K_w', 'P_c', 'mu_w', '1/dy' for node 2
    
    " Description of dictionary terms are as follows: "
    # K_w: permeability of porous domain [m^2]
    # P_c: capillary pressure = -2*gamma*cos(theta)/r [Pa]
    # mu_w: viscosity of liquid water [Pa-s]
        
    TDY1, TDY2 = TDY_vec[0], TDY_vec[1]
    eps_w1, eps_w2 = eps_vec[0], eps_vec[1]
    
    if i == 'channel':
        half = 2
        # r_p1 = p1['r_p'] *np.sqrt(p1['eps_g'][0] / p1['eps_go'][0])
        # r_p2 = p1['r_p'] *np.sqrt(p1['eps_g'][0] / p1['eps_go'][0])
        K_w = p.d['wt1']*p1['K_w'][0] + p.d['wt2']*p2['K_w'][0]
        
        if p2['ca_w'] < 90: J2 = J_hydrophilic(p2['s_w'])
        else: J2 = J_hydrophobic(p2['s_w'])
        
        P1_c = 0
        P2_c = P_cap(p2)[0]*J2[0]
        
    elif p.name == 'gdl_cl':
        half = 1
        # r_p1 = p1['r_p'] *np.sqrt(p1['eps_g'][-1] / p1['eps_go'][-1])
        # r_p2 = p2['r_p'] *np.sqrt(p2['eps_g'][0] / p2['eps_go'][0])
        K_w = p.d['wt1']*p1['K_w'][-1] + p.d['wt2']*p2['K_w'][0]
        
        if p1['ca_w'] < 90: J1 = J_hydrophilic(p1['s_w'])
        else: J1 = J_hydrophobic(p1['s_w'])
        
        if p2['ca_w'] < 90: J2 = J_hydrophilic(p2['s_w'])
        else: J2 = J_hydrophobic(p2['s_w'])
        
        P1_c = P_cap(p1)[-1]*J1[-1]
        P2_c = P_cap(p2)[0]*J2[0]
        
    else:
        half = 1
        # r_p1 = p1['r_p'] *np.sqrt(p1['eps_g'][i] / p1['eps_go'][i])
        # r_p2 = p2['r_p'] *np.sqrt(p2['eps_g'][i+1] / p2['eps_go'][i+1])
        K_w = p.d['wt1']*p1['K_w'][i] + p.d['wt2']*p2['K_w'][i+1]
        
        if p1['ca_w'] < 90: J1 = J_hydrophilic(p1['s_w'])
        else: J1 = J_hydrophobic(p1['s_w'])
        
        if p2['ca_w'] < 90: J2 = J_hydrophilic(p2['s_w'])
        else: J2 = J_hydrophobic(p2['s_w'])
        
        P1_c = P_cap(p1)[i]*J1[i]
        P2_c = P_cap(p2)[i+1]*J2[i+1]
                    
    # set state 1 properties:
    ca.gas.TDY = TDY1
    P1_g = ca.gas.P
    # P1_c = 2*p1['st_w']*np.cos(np.deg2rad(p1['ca_w'])) / r_p1
    P1_w = P1_g - P1_c
        
    # set state 2 properties:
    ca.gas.TDY = TDY2
    P2_g = ca.gas.P
    # P2_c = 2*p2['st_w']*np.cos(np.deg2rad(p2['ca_w'])) / r_p2
    P2_w = P2_g - P2_c    
    
    # water properties
    V_bar_w = ca.h2o_b.volume_mole
    
    # molar concentration at boundary
    C_w = (p.d['wt1']*eps_w1 + p.d['wt2']*eps_w2) / V_bar_w
    
    molar_flux = -C_w*K_w*(P2_w - P1_w)*half*p.d['1/dy'] / p1['mu_w']
    
    if eps_w1 <= 0: molar_flux = min(molar_flux, 0)
    if eps_w2 <= 0: molar_flux = max(molar_flux, 0)
        
    # Output returns molar flux [mol/m^2-s]
    return molar_flux

