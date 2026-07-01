import numpy as np
from nasa_data import NASA_DATA
from constants import R


def _get_coeffs(name, T):
    d = NASA_DATA[name]
    a = d["a_low"] if T < d["t_mid"] else d["a_high"]
    return a, d["molar_mass"]


def cp(name, T):
    a, M = _get_coeffs(name, T)
    return R * (a[0] + a[1]*T + a[2]*T**2 + a[3]*T**3 + a[4]*T**4) / M


def cv(name, T):
    a, M = _get_coeffs(name, T)
    cp_mole = R * (a[0] + a[1]*T + a[2]*T**2 + a[3]*T**3 + a[4]*T**4)
    return (cp_mole - R) / M


def enthalpy(name, T):
    a, M = _get_coeffs(name, T)
    h_mole = R * (a[0]*T + a[1]*T**2/2 + a[2]*T**3/3 + a[3]*T**4/4 + a[4]*T**5/5 + a[5])
    return h_mole / M


def mixture_cp(composition, T):
    cp_mole = 0.0
    for name, mole_frac in composition.items():
        a, _ = _get_coeffs(name, T)
        cp_mole += mole_frac * R * (a[0] + a[1]*T + a[2]*T**2 + a[3]*T**3 + a[4]*T**4)
    M_mix = sum(mole_frac * NASA_DATA[n]["molar_mass"] for n, mole_frac in composition.items())
    return cp_mole / M_mix


def mixture_cv(composition, T):
    cp_mole = 0.0
    for name, mole_frac in composition.items():
        a, _ = _get_coeffs(name, T)
        cp_mole += mole_frac * R * (a[0] + a[1]*T + a[2]*T**2 + a[3]*T**3 + a[4]*T**4)
    M_mix = sum(mole_frac * NASA_DATA[n]["molar_mass"] for n, mole_frac in composition.items())
    return (cp_mole - R) / M_mix


def mixture_enthalpy(composition, T):
    h_mole = 0.0
    for name, mole_frac in composition.items():
        a, _ = _get_coeffs(name, T)
        h_mole += mole_frac * R * (a[0]*T + a[1]*T**2/2 + a[2]*T**3/3 + a[3]*T**4/4 + a[4]*T**5/5 + a[5])
    M_mix = sum(mole_frac * NASA_DATA[n]["molar_mass"] for n, mole_frac in composition.items())
    return h_mole / M_mix
