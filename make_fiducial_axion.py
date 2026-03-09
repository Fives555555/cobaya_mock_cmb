#!/usr/bin/env python3
# Example script to create fiducial values for mock CMB likelihoods
from cobaya.model import get_model
import camb
import numpy as np
from camb.axion_utils import get_axion_phi_i

H0 = 67.32
h = H0 / 100

# omega_m = 0.3158
# baryon_fraction = 0.157

# MPc_in_sec**2 /Tpl**2, use to convert to 1/Mpc^2
# units = 1.4503545664380023e+113
# # eV to MPl
ev_conv = 1 / (2.435e27)
# ombh2 = omega_m * h**2 * baryon_fraction
# omch2 = omega_m * h**2 * (1 - baryon_fraction)
ombh2 = 0.022383
omch2 = 0.12011
# nnu = 3.046
As = 2.100e-09
ns = 0.96605
tau = 0.0543

accuracy = 2.5
# f_axion = 0.02
# mass = 1e-28
mH = 50
weighting_factor = 10.0

masses = [1e-28, 1e-27, 1e-26, 1e-25, 1e-24, 1e-23]
str_masses = [28, 27, 26, 25, 24, 23]
# TODO Change f_ax for each mass - current is just under upper 2 sigma constraints
# f_axs = [0.008, 0.02, 0.043, 0.58, 0.93, 0.99]
# All f_axs are 2% to match the forecast from the Simons Forecast paper
# f_axs = [0.02, 0.02, 0.02, 0.02, 0.02, 0.02]
# masses = [1e-24, 1e-23]
# str_masses = [24, 23]
# f_axs = [0.93, 0.99]
# masses = [1e-23]
# str_masses = [23]
# f_axs = [0.99]

# TODO 1e-23 case will not work, fix it

##########################################################

for j in range(len(masses)):
    mass = masses[j]
    f_axion = 0.02 #f_axs[j]
    accuracy_boost = accuracy
    axion_params_dict = None
    if f_axion > 0:
        print(f"\n0. Finding axion initial conditions for f_ax={f_axion}...")
        axion_params_dict = get_axion_phi_i(
            h=h,
            ombh2=ombh2,
            omch2_total=omch2,
            f_ax=f_axion,
            mass_ev=mass,
            mH=mH,
            use_PH=True,
            weighting_factor=weighting_factor,
            verbose=False,
            accuracy=accuracy_boost,
        )
        if axion_params_dict is None:
            raise RuntimeError("Failed to find axion initial conditions")

    omch2_cdm = (1 - f_axion) * omch2
    omch2_cdm = max(omch2_cdm, 1e-7)
    # camb:
    # from best fit with fixed massless neutrinos and nuisance-marginalized high-l
    fiducial_params = {
        # LambdaCDM parameters
        'H0': H0,
        # '100*theta_s': 1.041920539e+00,
        'ombh2': ombh2,
        # 'nnu': 3.046,  # three massless neutrinos
        'omch2': omch2_cdm,
        'As': As,
        # 'sigma8': 8.245006041e-01,
        'ns':  ns,
        'tau':  tau,
        
    }

    fiducial_params_extra = {
        'kmax': 1,
        'k_per_logint': 130,
        'AccuracyBoost':2.5,
        'lens_margin': 2050,
        'lAccuracyBoost': 1.2,
        'min_l_logl_sampling': 6000,
        'DoLateRadTruncation': False,
        # 'recombination': 'recfast',
        # 'halofit_version': 'mead',    
        'dark_energy_model': 'EarlyQuintessence',
        # 'AccuracyBoost': accuracy,
        'm': axion_params_dict['m'],
        'potential_type': 1,
        'theta_i': axion_params_dict['theta_i'],
        'mH': mH,
        'use_fluid_approximation': True,
        'use_PH': True,
        'weighting_factor': weighting_factor,
        'frac_lambda0': axion_params_dict['frac_lambda0'],
        'use_zc': False,
        'nonlinear': False,
        'lens_potential_accuracy': 0,
        'NonLinear': 'NonLinear_none'
    }


    fiducial_params_full = fiducial_params.copy()
    fiducial_params_full.update(fiducial_params_extra)

    info_fiducial = {
        'params': fiducial_params,
        'likelihood': {f'cobaya_mock_cmb.MockSOAx{str_masses[j]}': {'python_path': '.'},
                       f'cobaya_mock_cmb.MockPlanckLowLAx{str_masses[j]}': {'python_path': '.'}},
        'theory': {'camb': {"extra_args": fiducial_params_extra}}}

    model_fiducial = get_model(info_fiducial)

    model_fiducial.logposterior({})

    Cls = model_fiducial.provider.get_Cl(units="muK2")

    for likelihood in model_fiducial.likelihood.values():
        likelihood.create_fid_values(Cls, fiducial_params_full, override=True)
