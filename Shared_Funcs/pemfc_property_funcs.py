"""
Calculate Nafion transport properties as a function of morphology/state:
    
    The Nafion electrolyte used in the particle_shell_pemfc_* files has a
    complex microstructure and morphology that affects the value of important
    parameters that are used in the model, i.e. ionic conductivity and 
    effective oxygen diffusion coefficient.
    
    These parameters are functions of the microstructure of Nafion as well as 
    the state of the local environment, i.e. temp, RH, etc. The defined 
    functions below take information found in relavent literature and 
    approximate appropriate transport parameters based on the supplied user 
    inputs from the model.    
"""

""" Import helpful modules """
"-----------------------------------------------------------------------------"
import numpy as np

""" Nafion ionic conductivity model [S/m] """
"-----------------------------------------------------------------------------"
def RH_eq_func(rho_av,p,i):
    # Using excel, logarithm curves were fit to rho_n_av values at known RHs
    # to develop a function that takes in rho_n_av and returns an equivalent RH
    
    t_yamls = np.array([5,7,12,18,20,42])
    t_i = np.argmin(abs(t_yamls - p['t_naf'][i]*1e9))
    
    C0 = np.array([0.2105, 0.1987, 0.2318, 0.2270, 0.2290, 0.2281])
    C1 = np.array([0.1938, 0.1659, 0.2361, 0.2274, 0.2310, 0.2290])
    
    RH_eq = C0[t_i]*np.log(rho_av) - C1[t_i]            
    
    return RH_eq

def sig_io_func(p, ca, sv):
    # The method input gives control over how the Nafion conductivity is 
    # calculated. Options are 'lam' for laminar in which an interpolation is
    # done using data from [1], 'bulk' for treating the thin Nafion shells the 
    # as a bulk-like material using NR results from [5], and 'mix' which uses a
    # weighted parallel mixutre of 'lam' and 'bulk' based on how much Pt vs C 
    # exists at current conditions. This is because it is speculated that Pt 
    # may have lamellae although C may not. 'sun' was also added to the 
    # agglomerate model options which takes constant values used in [2].
    
    # Inputs: Temperature [K], Nafion shell thickness [m], rel. humiditiy [%],
    #         Pt coverage [%], p['eps/tau2_n'] [-] and p['p_eff_SAnaf'] [-], 
    #         and calculation method [-]
    
    """ Lamellae Method """
    # Data below is taken from "Proton Transport in Supported Nafion Nanothin
    # Films by Electrochemical Impedence Spectroscopy" by Paul, MacCreery, and
    # Karan in their Supporting Information Document [1]. The data was given in 
    # mS/cm and converted to S/m for the model calling this function.
    
    # indecies: temperature [C], Nafion shell thickness [nm], RH [%] 
    sig_data = np.zeros([5,5,5]) 
    temp_vals = np.array([25,30,40,50,60])
    thick_vals = np.array([4,10,55,160,300])
    RH_vals = np.array([20,40,60,80,95])
        
    # v_w = np.zeros([p['Ny'],p['Nr']])
    # for i in range(p['Ny']):
        
    #     ih_n = ca.naf_b[i].species_index('H(Naf)')
    #     ih2o_n = ca.naf_b[i].species_index('H2O(Naf)')
        
    #     for j in range(p['Nr']):
    #         ca.naf_b[i].Y = sv[ca.ptr['rho_naf_k'] +i*p['nxt_y'] +j*p['nxt_r']]
            
    #         v_k = ca.naf_b[i].X*ca.naf_b[i].partial_molar_volumes
    #         v_w[i,j] =  v_k[ih2o_n] / sum(v_k)
    
    # v_w_a = np.sum(p['Vf_shl']*v_w,axis=1)
    # lamb_n = np.clip((v_w_a / (1 - v_w_a) *983/1980 *1100/18.02), 0., 22.)
    
    rho_naf_w = np.zeros([p['Ny'],p['Nr']])
    for i in range(p['Ny']):
        
        ih2o_n = ca.naf_b[i].species_index('H2O(Naf)')
        
        for j in range(p['Nr']):
            ind = ca.ptr['rho_naf_k'] +i*p['nxt_y'] +j*p['nxt_r']
            rho_naf_w[i,j] = sv[ind][ih2o_n]    
    
    rho_naf_av = np.sum(p['Vf_shl']*rho_naf_w,axis=1)
    
    RH, RH_C = np.zeros(p['Ny']), np.zeros(p['Ny'])
    for i in range(p['Ny']):
        av = rho_naf_av[i]
        if av > 0: 
            RH[i] = RH_eq_func(av,p,i)*100
            RH_C[i] = RH_eq_func(av/2,p,i)*100
        else: 
            RH[i] = min(RH_vals)
            RH_C[i] = min(RH_vals)
    
    "Data for 25C as thickness[nm] for rows and RH[%] for columns"
    sig_data[0,:,:] = np.array([[0.0002,0.0206,0.4138,4.9101,21.888],   # t 4nm
                                [0.0002,0.0199,0.4073,5.1758,23.9213],  # t 10nm
                                [0.0002,0.0269,0.5448,5.3493,22.753],   # t 55nm
                                [0.3362,3.2505,8.3065,27.0725,54.0428], # t 160nm
                                [1.5591,8.8389,19.6728,None,None]])     # t 300nm
    
    "Data for 30C as thickness[nm] for rows and RH[%] for columns"
    sig_data[1,:,:] = np.array([[0.0001,0.012,0.278,3.432,21.481],   # t 4nm
                                [0.0003,0.018,0.339,3.895,22.062],   # t 10nm
                                [0.0004,0.028,0.550,4.296,20.185],   # t 55nm
                                [0.0016,0.081,1.120,9.244,34.810],   # t 160nm
                                [0.0071,0.359,2.797,10.978,43.913]]) # t 300nm
    
    "Data for 40C as thickness[nm] for rows and RH[%] for columns"
    sig_data[2,:,:] = np.array([[0.0003,0.029,0.585,6.164,30.321],   # t 4nm
                                [0.0009,0.034,0.625,5.374,48.799],   # t 10nm
                                [0.0011,0.065,0.931,6.909,40.439],   # t 55nm
                                [0.0032,0.152,1.770,14.162,68.326],  # t 160nm
                                [0.0140,0.605,4.939,17.083,68.334]]) # t 300nm
    
    "Data for 50C as thickness[nm] for rows and RH[%] for columns"
    sig_data[3,:,:] = np.array([[0.001,0.062,1.087,8.335,37.686],   # t 4nm
                                [0.002,0.077,1.031,8.127,57.339],   # t 10nm
                                [0.002,0.121,1.603,9.149,48.934],   # t 55nm
                                [0.007,0.247,2.704,19.221,72.006],  # t 160nm
                                [0.031,1.076,7.185,20.981,83.923]]) # t 300nm
    
    "Data for 60C as thickness[nm] for rows and RH[%] for columns"
    sig_data[4,:,:] = np.array([[0.003,0.14,1.51,11.16,55.18],   # t 4nm
                                [0.003,0.17,1.72,13.67,62.39],   # t 10nm
                                [0.007,0.24,2.29,16.60,63.20],   # t 55nm
                                [0.015,0.45,4.31,26.63,93.33],   # t 160nm
                                [0.009,0.44,3.43,26.73,100.60]]) # t 300nm
    
    "Create interpolation function for relavent ranges"
    from scipy.interpolate import RegularGridInterpolator
    sig_io_int = RegularGridInterpolator((temp_vals,thick_vals,RH_vals),sig_data)
    
    "Call interpolation function for model specified paramaters"
    # Multiplication by 0.1 is unit conversion from mS/cm to S/m. Runner file
    # stores T and t_naf in [K] and [m] so are also converted inside the 
    # interpolation function to the same units as original data [C] and [nm].
    RH = np.clip(RH, min(RH_vals), max(RH_vals))  
    RH_C = np.clip(RH_C, min(RH_vals), max(RH_vals))

    pts = np.zeros([p['Ny'],3])
    for i in range(p['Ny']):
        pts[i,:] = [p['T']-273, p['t_naf'][i]*1e9, RH[i]]
        
    sig_io_lam = sig_io_int(pts) *0.1
    
    """ Bulk Method """
    # This method assumes that the thin shell of Nafion is treated the same as
    # the bulk material. Lambda is calculated using an empirical relationship. 
    # Then the sig_io formula from [5] for a bulk membrane is used and scaled 
    # by the scaling factor, also from [5].
    
    # The loop below assumes RH is not RH_eq and instead is the actual local 
    # gas-phase RH. 
    
    if p['sig_method'] == 'lit':
        for i in range(p['Ny']):
            ih2o_g = ca.gas.species_index('H2O')
            
            rho_gas_k = sv[ca.ptr['rho_gas_k'] +i*p['nxt_y']]
            ca.gas.TDY = p['T'], sum(rho_gas_k), rho_gas_k
            
            RH[i] = ca.gas.X[ih2o_g]*ca.gas.P / 19946 *100
    
    lamb_n = 0.3 + 10.8*(RH/100) - 16*(RH/100)**2 + 14.1*(RH/100)**3
    sig_io_lit = (0.5139*lamb_n - 0.326)*np.exp(1268*(1/303 - 1/p['T']))
    sig_io_bulk = sig_io_lit *0.672
        
    """ Mix Method """
    # Using a parallel resistor network to weight the conductivity through 
    # lamellae and that through bulk-like material is performed with respect to
    # the amount of Pt and C areas respectively.
    sig_io_mix = 1 / (p['p_Pt']/100 /sig_io_lam +(1-p['p_Pt']/100) /sig_io_bulk)
    
    " Set conductivity depending on method "
    # Based on the method, return the appropriate conductivity.
    if p['sig_method'] == 'lam': sig_io = sig_io_lam 
    elif p['sig_method'] == 'bulk': sig_io = sig_io_bulk 
    elif p['sig_method'] == 'mix': sig_io = sig_io_mix 
    elif p['sig_method'] == 'lit': sig_io = sig_io_lit
    
    # Output returns ionic conductivity [S/m]
    return sig_io

""" Effective O2 diffusion coeff. in Nafion model [m^2/s] """
"-----------------------------------------------------------------------------"
def D_naf_k_func(p):
    # The method input gives control over how the Nafion conductivity is 
    # calculated. Options are 'lam' for laminar in which artificial lamellae
    # are generated from the shell thickness and used to created a series 
    # resistance network to approximate D_O2. Lamellae thicknesses, water
    # volume fractions, and scalings are taken from [5]. Additionally, 'bulk' 
    # can be used in order to treat the thin Nafion shell the same as bulk 
    # material using fits from [4] and water scaling from [5]. 'sun' was also 
    # added to the agglomerate model options which takes constant values used 
    # in [2].
    
    # Inputs: Temperature [K], Carbon radius [m], Nafion shell thickness [m],
    #         Pt coverage [%], p['p_eff_SAnaf'] [-], calculation method [-],
    #         and model type [-]
    
    t_naf = np.mean(p['t_naf'])
    
    """ Lamellae Method """ 
    # This method assumes that lamellae exist in the thin Nafion shells found
    # in the CL. A series resistor network in spherical coordinates is used to
    # approximate the effective O2 diffusion coefficient in this case. Starting
    # with a bulk O2 value, each lamellae is scaled by V_w_i and f_i taken from
    # values used in conductivity fits in [5].
    
    "Data from lamellae thicknesses [A], V_w [-], and f [-] taken from [5]"
    # Scaling factors are in order of f_lam_1, f_lam_2, bulk_like... where each
    # lamellae has a scaling factor that is linear between f_lam_1 and f_lam_2
    # as a function of its layer number.
    f = np.array([0.33456564, 0.7488917, 0.671766])
    lam_num = np.array([4, 5, 6, 5])    
    
    t_i = np.array([[8.612, 24.118, 10.505, 15.281, 0, 0, 0],
                    [7.5947, 20.907, 15.228, 18.3322, 26.467, 0, 0],
                    [10.866, 16.236, 16.600, 16.202, 21.114, 11.830, 417.886],
                    [9.481, 22.070, 16.776, 19.297, 18.831, 1477.106, 56.649]])
    
    V_w_i = np.array([[0.659, 0.108, 0.52138368, 0.05217589, 0., 0., 0.],
                      [0.659, 0.108, 0.33375568, 0.07767769, 0.31938309, 0., 0.],
                      [0.659, 0.108, 0.27107169, 0.08195979, 0.19785204, 0.13540123, 0.17599638],
                      [0.81347793, 0.08000953, 0.57073381, 0.19147129, 0.33532437, 0.25710187, 0.20889646]])
    
    "Use arrhenius relationship from Sethuraman et al. [4] to get ref D_O2"
    D_eff_ref = 17.45e-6*np.exp(-1514 /p['T']) /100**2 # 100**2 [cm^2/s -> m^2/s]
    
    "Scale each lamellae's D_O2 by V_w_i and f_i"
    D_O2_i = D_eff_ref *(V_w_i /0.367)
    
    for i in range(t_i.shape[0]):
        f_lam = np.linspace(f[0], f[1], lam_num[i])
        
        if t_i[i, lam_num[i]] != 0:
            f_bulk = np.ones(t_i.shape[1] - lam_num[i])*f[2]
        elif t_i[i, lam_num[i]] == 0:
            f_bulk = np.zeros(t_i.shape[1] - lam_num[i])
            
        f_i = np.hstack([f_lam, f_bulk])
        D_O2_i[i] = D_eff_ref *(V_w_i[i,:] /0.367) *f_i
        
    "Build series resistor network to get total effective D_O2"
    R_i = np.zeros_like(t_i)
    
    for i in range(t_i.shape[0]):        
        for j in range(np.count_nonzero(t_i[i,:])):
            R_i[i,j] = 1 / D_O2_i[i,j]
            
    R_avg = np.sum(t_i*R_i, axis=1) / np.sum(t_i, axis=1)
    D_O2_lam_vec = 1 / R_avg
    
    "Interpolate between film thicknesses to get approximate D_O2"
    D_eff_naf_lam = np.interp(t_naf, np.sum(t_i, axis=1) /10 *1e-9, D_O2_lam_vec)
    
    """ Bulk Method """
    # Relationships below use temperature to approximate an oxygen diffusion
    # coefficient in bulk Nafion from [4]. This value is then normalized by
    # the water volume from their samples and rescaled by the assumed water
    # volume in the shell based on t_naf and data taken from [5].
    
    "Data from [5] in [A] for 12, 18, 20, 42, 60, 103, 120, 154nm t_bulk"
    # Bulk-like thicknesses from [5] for all samples. First two were lamellae 
    # only. Also, water volme fractions for the bulk-like layers. Although 
    # original data was in [A], convertion to [m] is done below.
    tb_i = np.array([44.316, 123.3766, 143.297, 417.886, 720.537, 1257.711 +35.4610,
                     1477.106 +56.649, 1992.758]) / 10 *1e-9
    V_wb_i = np.array([0.122, 0.152, 0.147, 0.176, 0.246, 0.246, 0.258, 0.255])
    
    "Use arrhenius relationship from Sethuraman et al. [4] to get ref D_O2"
    D_eff_ref = 17.45e-6*np.exp(-1514 /p['T']) /100**2 # 100**2 [cm^2/s -> m^2/s]    
    
    "Find interpolation of V_wb_i for specified t_naf"
    if t_naf < tb_i[0]:
        V_w_interp = V_wb_i[0]
    else:
        V_w_interp = np.interp(t_naf, tb_i, V_wb_i)
    
    "Normalize V_w with 0.367 - from lambda=18 in [4]"
    D_eff_naf_bulk = D_eff_ref *(V_w_interp /0.367)
    
    """ Mix Method """
    # Using a parallel resistor network to weight the diffusion coeff. through 
    # lamellae and that through bulk-like material is performed with respect to
    # the amount of Pt and C areas respectively.
    D_eff_naf_mix = D_eff_naf_lam
    
    " Set diffusion coefficient depending on method "
    if p['D_O2_method'] == 'lam':
        D_eff_naf = D_eff_naf_lam
    elif p['D_O2_method'] == 'bulk':
        D_eff_naf = D_eff_naf_bulk
    elif p['D_O2_method'] == 'mix':
        D_eff_naf = D_eff_naf_mix
    elif p['D_O2_method'] == 'sun':
        D_eff_naf = 8.45e-10 # from [2]
    
    # Output returns effective O2 diffusion coefficient [m^2/s]
    return D_eff_naf

"""    
Offload geometric calculations for reaction areas to clean up the code:
    
    Using user input for carbon and Pt radii as well as a Pt-loading (common
    in PEMFC fabrication), the reaction areas are estimated. These calculations
    were combined into this set of functions in order to provide a shorter code
    for the runner and dsvdt functions therefore making them more managable to
    edit and debug.
"""

""" Caclulate the reaction areas from Pt-loading and radii information """
# The following reaction areas function is defined for the core-shell micro-
# structure. In the function name, "cs" is used to stand for "core-shell."
"-----------------------------------------------------------------------------"
def rxn_areas_cs(w_Pt,t_cl,eps_gas,t_naf,r_c,r_Pt,rho_Pt,theta,Ny_cl,a_n):
    # Units for inputs are:
    # w_Pt [mg/cm^2], t_cl [m], eps_gas [-], t_naf [m], r_c [m], r_Pt [m], 
    # rho_Pt [kg/m^3], theta [degrees], Ny_cl [-], a_n [-]
    
    dy_cl = t_cl / Ny_cl
    
    "Max t_naf for Pt distribution calc"
    t_naf = max(t_naf)
    
    "Find the mass of Pt per agglomerate"
    w_Pt = w_Pt *0.01 # convert [mg/cm^2] --> [kg/m^2]
    V_part = 4/3*np.pi*(r_c + t_naf)**3
    m_Pt_part = a_n.d['wt']*w_Pt / dy_cl / (1 - eps_gas) *V_part
    
    "Distribute Pt mass to half sphere bulbs"
    V_Pt_part = m_Pt_part / rho_Pt
    V_Pt_bulb = 2/3*np.pi*r_Pt**3 # volume of Pt 1/2 sphere sitting on C surf
    N_bulbs = V_Pt_part / V_Pt_bulb
    
    "Using r_Pt and assuming semi-spheres, find Pt surface area"
    SA_c_part = 4*np.pi*r_c**2
    V_c_part = 4/3*np.pi*r_c**3
    
    SA_Pt_part = N_bulbs*2*np.pi*r_Pt**2
    SA_dl_part = SA_c_part - N_bulbs*np.pi*r_Pt**2 + SA_Pt_part
    SA_naf_int = 4*np.pi*(r_c + t_naf)**2
    V_naf_part = V_part - V_c_part - V_Pt_part
    
    "Determine method for Nafion SA based on theta"
    p_Pt = SA_Pt_part / SA_dl_part
    SA_naf_edit = 4*np.pi*r_c**2*p_Pt*(1+np.tan(np.deg2rad(theta)))**2
    
    SA_pv_naf = np.zeros_like(SA_naf_edit)
    p_eff_SAnaf = np.zeros_like(SA_naf_edit)
    
    for i in range(SA_pv_naf.size):
        
        if SA_naf_int < SA_naf_edit[i]:
            SA_pv_naf[i] = (1 - eps_gas)*SA_naf_int / V_part
            p_eff_SAnaf[i] = 1.0
        else:
            SA_pv_naf[i] = (1 - eps_gas)*SA_naf_edit[i] / V_part
            p_eff_SAnaf[i] = SA_naf_edit[i] / SA_naf_int
        
    "Use SA and V factors to find geometric parameters to return"
    geom_out ={}
    geom_out['SApv_naf'] = SA_pv_naf # Nafion SA per V
    geom_out['SApv_pt'] = (1 - eps_gas)*SA_Pt_part / V_part # Pt SA per V
    geom_out['SApv_dl'] = (1 - eps_gas)*SA_dl_part / V_part # Double layer SA per V
    geom_out['eps_n_cl'] = (1 - eps_gas)*V_naf_part / V_part # Nafion V fraction
    geom_out['p_eff_SAnaf'] = p_eff_SAnaf # % effective SA of Nafion between shells
    geom_out['V_c_part'] = V_c_part
    geom_out['V_Pt_part'] = V_Pt_part
    geom_out['SA_dl_part'] = SA_dl_part
    geom_out['N_pt_part'] = N_bulbs
  
    # Units for outputs are:
    # SA_pv_* [1/m], eps_n_cl [-], p_eff_SAnaf [-]
    return geom_out
    