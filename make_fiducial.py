#!/usr/bin/env python3
# Example script to create fiducial values for mock CMB likelihoods
from cobaya.model import get_model

# # class
# # from best fit with fixed massless neutrinos and nuisance-marginalized high-l
# fiducial_params = {
#     # LambdaCDM parameters
#     'H0': 68.1461,#67.66,
#     # '100*theta_s': 1.041920539e+00,
#     'omega_b': 2.241606328e-02,#0.02242,
#     'N_ur': 3.046,  # three massless neutrinos
#     'omega_cdm': 1.194462584e-01,#0.11933,
#     'A_s': 2.109371924e-09,#2.105e-09,
#     # 'sigma8': 8.245006041e-01,
#     'n_s': 9.660360599e-01,#0.9665,
#     'tau_reio': 5.142494234e-02#0.0561
# }

# fiducial_params_extra = {
#     'recombination': 'recfast',
#     'non linear': 'halofit'
# }
# # class


#TODO PUT THE AXION PARAMS CALCULATION IN

# H0 = 67.32
# h = H0 / 100

# omega_m = 0.3158
# baryon_fraction = 0.157

# # MPc_in_sec**2 /Tpl**2, use to convert to 1/Mpc^2
# units = 1.4503545664380023e+113
# # eV to MPl
# ev_conv = 1 / (2.435e27)
# # ombh2 = omega_m * h**2 * baryon_fraction
# # omch2 = omega_m * h**2 * (1 - baryon_fraction)
# ombh2 =  0.022383
# omch2 =  0.12011

# lmax = 3500
# accuracy=2.5
# f_axion = 1
# mass = 1e-28
# mH = 150
# weighting_factor = 10.0

# # pars = camb.set_params(H0=67.32, ombh2=0.022383, omch2=0.12011, mnu=0.06, omk=0, tau=0.0543, As=np.exp(3.0448)/1e10, ns=0.96605, halofit_version='mead', lmax=3500)

# # To match the SO mock fiducial:
# # pars = camb.set_params(H0=68.1461, ombh2=0.02241606328, nnu=3.046, omch2=0.1194462584, As=2.109371924e-09, ns=0.9660360599, tau=0.05142494234, halofit_version='mead', lmax=3500)

# # To match Planck 2018 (comaprison to fiducial used currently):
# # pars = camb.set_params(H0=67.66, ombh2=0.02242, nnu=3.046, omch2=0.11933, As=2.105e-09, ns=0.9665, tau=0.0561, halofit_version='mead', lmax=3500, AccuracyBoost=2.5)



# #####################################
# # To include an axion in the fiducial model (using to compare to that made with the make_fiducial.py file)
# # determine correct background NO AXIONS
# pars = camb.CAMBparams()
# pars = camb.set_params(H0=H0, ombh2=ombh2, omch2=omch2, lmax=lmax, AccuracyBoost=accuracy)

# results = camb.get_background(pars)
# frac_lambda0 = results.grhov / (results.grhov + f_axion * results.grhoc)

# de = camb.CAMBparams.make_class_named(
#     camb.dark_energy.EarlyQuintessence, 
#     camb.dark_energy.DarkEnergyModel)

# # Update omch2
# new_omch2 = (1 - f_axion) * omch2
# pars = camb.set_params(H0=H0, ombh2=ombh2, omch2=new_omch2, lmax=lmax, AccuracyBoost=accuracy)

# omaxh2 = f_axion * omch2              
# phi_i = 1e-5
# factor = 1.5


# def get_rho_ax_crit(results, z):
#     a = 1 / (1 + z)
#     rho, _ = results.get_dark_energy_rho_w(a)
#     return (rho - results.frac_lambda0*results.grhov) / results.grhocrit * a**3

# # determine correct phi INCLUDE AXIONS

# for i in range(100):
#     phi_i = phi_i * factor
#     de.set_params(m=mass*ev_conv, potential_type=1, theta_i=phi_i, 
#     mH=mH, use_fluid_approximation=True, use_PH=True, 
#         weighting_factor=weighting_factor, frac_lambda0=frac_lambda0, use_zc=False)
#     pars.DarkEnergy = de
#     results = camb.get_background(pars)        
#     if get_rho_ax_crit(results, 0.0) * (pars.H0/100)**2 > omaxh2:
#         break
#     if i == 99:
#         print("Failed to find phi_i (stage 1)")
# phi_i_lower = phi_i / factor
# phi_i_upper = phi_i
# for i in range(100):
#     phi_i = (phi_i_lower + phi_i_upper) / 2
#     de.set_params(m=mass*ev_conv, potential_type=1, theta_i=phi_i, 
#         mH=mH, use_fluid_approximation=True, use_PH=True, 
#             weighting_factor=weighting_factor, frac_lambda0=frac_lambda0, use_zc=False)
#     pars.DarkEnergy = de
#     results = camb.get_background(pars)      
#     if get_rho_ax_crit(results, 0.0) * (pars.H0/100)**2 > omaxh2:
#         phi_i_upper = phi_i
#     else:
#         phi_i_lower = phi_i
#     if np.abs(get_rho_ax_crit(results, 0.0) * (pars.H0/100)**2 - omaxh2) < 1e-4:
#         # print(f'Phi_i = {phi_i:.3e}, omaxh2 = {omaxh2:.3e}, f_axion = {f_axion:.3e}')
#         break
#     if i == 99:
#         raise print("Failed to find phi_i (stage 2)")

# de.set_params(m=mass*ev_conv, potential_type=1, theta_i=phi_i, 
#     mH=mH, use_fluid_approximation=True, use_PH=True, 
#         weighting_factor=weighting_factor, frac_lambda0=frac_lambda0, use_zc=False)
# pars.DarkEnergy = de   
# ###############################################











# camb:
# from best fit with fixed massless neutrinos and nuisance-marginalized high-l
fiducial_params = {
    # LambdaCDM parameters
    'H0': 67.32,
    # '100*theta_s': 1.041920539e+00,
    'ombh2': 0.022383,
    'nnu': 3.046,  # three massless neutrinos
    'omch2': 0.12011,
    'As': 2.100e-09,
    # 'sigma8': 8.245006041e-01,
    'ns':  0.96605,
    'tau':  0.0543,
    
}

# CHANGE THE ELL MAX IN THE BASELINE.YAML TO MATCH THE NOISE CURVES


fiducial_params_extra = {
    'AccuracyBoost':2.5
    # 'recombination': 'recfast',
    # 'halofit_version': 'mead'
}
# camb


fiducial_params_full = fiducial_params.copy()
fiducial_params_full.update(fiducial_params_extra)

info_fiducial = {
    'params': fiducial_params,
    'likelihood': {'cobaya_mock_cmb.MockSO': {'python_path': '.'},
                   'cobaya_mock_cmb.MockSOBaseline': {'python_path': '.'},
                   'cobaya_mock_cmb.MockSOGoal': {'python_path': '.'},
                   'cobaya_mock_cmb.MockCMBS4': {'python_path': '.'},
                   'cobaya_mock_cmb.MockCMBS4sens0': {'python_path': '.'},
                   'cobaya_mock_cmb.MockPlanck': {'python_path': '.'},
                   'cobaya_mock_cmb.MockPlanckHighL': {'python_path': '.'},
                   'cobaya_mock_cmb.MockPlanckLowL': {'python_path': '.'}},
    #'theory': {'class': {"extra_args": fiducial_params_extra}}}
    'theory': {'camb': {"extra_args": fiducial_params_extra}}}

model_fiducial = get_model(info_fiducial)

model_fiducial.logposterior({})

Cls = model_fiducial.provider.get_Cl(units="muK2")

for likelihood in model_fiducial.likelihood.values():
    likelihood.create_fid_values(Cls, fiducial_params_full, override=True)
