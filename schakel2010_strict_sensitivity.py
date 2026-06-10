#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Standalone Schakel & Smeulders (2010) reproduction.

The script intentionally keeps the Schakel (2010) reproduction separate from
the finite-offset waveform model. It implements the Table I parameters,
Appendix A dynamic coefficients, Eq. (24)-(39) wave properties, Appendix B
boundary-value system, Fig. 2 potential-coefficient angle scan, and Eq. (47)-(51)
vertical energy-flux coefficients.
"""

from __future__ import annotations

import argparse
import math
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Dict, Mapping, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@dataclass
class SchakelConfig:
    """Schakel & Smeulders (2010) Table I parameters."""

    K_b: float = 5.8e9
    G: float = 3.4e9
    K_s: float = 40.0e9
    K_f: float = 2.22e9
    eta: float = 1.0e-3
    rho_f: float = 1000.0
    rho_s: float = 2760.0
    alpha_inf: float = 2.3
    phi: float = 0.24
    k0_m2: float = 0.390e-12
    temperature: float = 295.0
    pH: float = 7.0
    eps_f: float = 80.0
    eps_s: float = 4.0
    z1: int = 1
    z2: int = -1
    b1: float = 3.246e11
    b2: float = 4.931e11
    M_similarity: float = 1.0
    C_molL: float = 1.0e-3

    K_fl: float = 2.22e9
    rho_fl: float = 1000.0
    eps_fl: float = 80.0
    sigma_fl: float = 5.0e-3

    eps0: float = 8.854187817e-12
    mu0: float = 4.0 * math.pi * 1e-7
    e_charge: float = 1.602176634e-19
    k_B: float = 1.380649e-23
    N_A: float = 6.02214076e23
    eps_complex: float = 1.0e-30


TABLE_II_REFERENCE = {
    10.0: {
        "TE_Pf_Pf": 4.6700e-1,
        "TE_Ps_Ps": 1.8083e-3,
        "TE_TM_TM": -2.6500e-10,
        "TE_SV_SV": 1.3084e-1,
        "TE_Pf_Ps": -2.0028e-10,
        "TE_Pf_TM": -2.7470e-12,
        "TE_Pf_SV": 7.4887e-11,
        "TE_Ps_TM": 6.7059e-10,
        "TE_Ps_SV": -1.0420e-8,
        "TE_TM_SV": 4.3970e-14,
        "RE_EE": 7.4497e-11,
        "RE_Pr_Pr": -4.0035e-1,
    },
    500_000.0: {
        "TE_Pf_Pf": 4.7334e-1,
        "TE_Ps_Ps": 1.3026e-1,
        "TE_TM_TM": -4.9895e-7,
        "TE_SV_SV": 7.0295e-2,
        "TE_Pf_Ps": -8.0868e-7,
        "TE_Pf_TM": -2.5059e-9,
        "TE_Pf_SV": 3.0585e-7,
        "TE_Ps_TM": 5.2570e-7,
        "TE_Ps_SV": -3.7939e-5,
        "TE_TM_SV": 1.0814e-7,
        "RE_EE": 1.1871e-7,
        "RE_Pr_Pr": -3.2615e-1,
    },
}

TABLE_III_REFERENCE = {
    30.0: {
        "TE_Pf_Pf": 4.7389e-1,
        "TE_Ps_Ps": 1.2478e-1,
        "TE_TM_TM": -5.0880e-7,
        "TE_SV_SV": 7.2674e-2,
        "TE_Pf_Ps": -2.6036e-6,
        "TE_Pf_TM": -3.3667e-9,
        "TE_Pf_SV": 9.8734e-7,
        "TE_Ps_TM": 7.1950e-7,
        "TE_Ps_SV": -1.2373e-4,
        "TE_TM_SV": 1.3113e-7,
        "RE_EE": 1.3503e-7,
        "RE_Pr_Pr": -3.2877e-1,
    },
    45.0: {
        "TE_Pf_Pf": 1.6454e-7,
        "TE_Ps_Ps": 2.7971e-1,
        "TE_TM_TM": -1.6581e-6,
        "TE_SV_SV": 6.6128e-1,
        "TE_Pf_Ps": 3.2114e-5,
        "TE_Pf_TM": 4.4496e-9,
        "TE_Pf_SV": 5.2354e-5,
        "TE_Ps_TM": 2.4385e-6,
        "TE_Ps_SV": 1.9413e-3,
        "TE_TM_SV": 8.4593e-7,
        "RE_EE": 5.4312e-7,
        "RE_Pr_Pr": -5.6982e-2,
    },
}


def complex_sqrt_branch(x: complex | np.ndarray) -> complex | np.ndarray:
    """Square-root branch for exp(i omega t) and exp(-i k z) decay."""
    y = np.sqrt(x + 0j)
    if np.isscalar(y):
        if np.imag(y) > 0 or (abs(np.imag(y)) < 1e-18 and np.real(y) < 0):
            y = -y
        return y
    mask = (np.imag(y) > 0) | ((np.abs(np.imag(y)) < 1e-18) & (np.real(y) < 0))
    y[mask] = -y[mask]
    return y


def h_concentration_for_ph(pH: float) -> float:
    """Return H+ concentration in mol/L for a requested pH."""
    return 10.0 ** (-float(pH))


def schakel_reference_config() -> Tuple[SchakelConfig, Dict[str, float]]:
    cfg = SchakelConfig()
    baseline = {
        "phi": cfg.phi,
        "k0_m2": cfg.k0_m2,
        "alpha_inf": cfg.alpha_inf,
        "cH_molL": h_concentration_for_ph(cfg.pH),
        "C_molL": cfg.C_molL,
        "omega": 1.0e6,
        "theta_deg": 45.0,
    }
    return cfg, baseline


def electrochemistry(C_molL: float, cH_molL: float, cfg: SchakelConfig) -> Dict[str, float]:
    C = float(C_molL)
    pH = -math.log10(max(float(cH_molL), cfg.eps_complex))
    zeta = (0.010 + 0.025 * math.log10(C)) * (pH - 2.0) / 5.0
    return {"C_molL": C, "cH_molL": cH_molL, "pH": pH, "zeta": zeta}


def dynamic_coefficients(
    phi: float,
    k0: float,
    alpha_inf: float,
    cH_molL: float,
    omega: float,
    cfg: SchakelConfig,
    C_override_molL: float | None = None,
) -> Dict[str, complex]:
    """Schakel Appendix A dynamic permeability, coupling, and conductivity."""
    C = cfg.C_molL if C_override_molL is None else float(C_override_molL)
    ec = electrochemistry(C, cH_molL, cfg)
    zeta = ec["zeta"]
    z_vals = np.array([cfg.z1, cfg.z2], dtype=float)
    b_vals = np.array([cfg.b1, cfg.b2], dtype=float)
    n_each = C * 1000.0 * cfg.N_A
    n_vals = np.array([n_each, n_each], dtype=float)

    omega_t = phi * cfg.eta / (alpha_inf * k0 * cfg.rho_f)
    Lambda = math.sqrt(max(8.0 * alpha_inf * k0 / (phi * cfg.M_similarity), cfg.eps_complex))
    k_dyn = k0 / (
        np.sqrt(1.0 + 1j * (omega / omega_t) * cfg.M_similarity / 2.0)
        + 1j * (omega / omega_t)
    )

    denom_debye = np.sum((cfg.e_charge * z_vals) ** 2 * n_vals) / (
        cfg.eps0 * cfg.eps_f * cfg.k_B * cfg.temperature
    )
    debye_d = math.sqrt(1.0 / max(float(np.real(denom_debye)), cfg.eps_complex))

    one_minus = 1.0 - 2.0 * debye_d / Lambda
    L = -(phi / alpha_inf) * (cfg.eps0 * cfg.eps_f * zeta / cfg.eta) * one_minus
    L *= (
        1.0
        + 2j
        * omega
        / (cfg.M_similarity * omega_t)
        * one_minus**2
        * (1.0 + debye_d * np.sqrt(1j * omega * cfg.rho_f / cfg.eta)) ** 2
    ) ** (-0.5)

    sigma_f = float(np.sum((cfg.e_charge * z_vals) ** 2 * b_vals * n_vals))
    if abs(zeta) < 1e-12 or sigma_f <= 0.0:
        C_em = 0.0
        C_os = 0.0 + 0j
        P_os = 0.0
    else:
        exp_terms = np.exp(-(cfg.e_charge * z_vals * zeta) / (2.0 * cfg.k_B * cfg.temperature)) - 1.0
        C_em = float(2.0 * debye_d * np.sum((cfg.e_charge * z_vals) ** 2 * b_vals * n_vals * exp_terms))
        P_os = float(
            (8.0 * cfg.k_B * cfg.temperature * debye_d**2)
            / (cfg.eps0 * cfg.eps_f * zeta**2)
            * np.sum(n_vals * exp_terms)
        )
        C_os = (
            0.0 + 0j
            if abs(P_os) < 1e-30
            else ((cfg.eps0 * cfg.eps_f) ** 2 * zeta**2 / (2.0 * debye_d * cfg.eta))
            * P_os
            * (1.0 + (2.0 / P_os) * debye_d * np.sqrt(1j * omega * cfg.rho_f / cfg.eta)) ** (-1)
        )

    sigma = (phi * sigma_f / alpha_inf) * (
        1.0 + 2.0 * (C_em + C_os) / max(sigma_f * Lambda, cfg.eps_complex)
    )
    eps_bulk = cfg.eps0 * (phi * (cfg.eps_f - cfg.eps_s) / alpha_inf + cfg.eps_s)
    eps_bar = eps_bulk - 1j * sigma / omega + 1j * cfg.eta * L**2 / (omega * k_dyn)

    return {
        **ec,
        "omega_t": omega_t,
        "Lambda": Lambda,
        "debye_d": debye_d,
        "k_dyn": k_dyn,
        "L": L,
        "sigma_f": sigma_f,
        "C_em": C_em,
        "C_os": C_os,
        "sigma": sigma,
        "eps_bulk": eps_bulk,
        "eps_bar": eps_bar,
    }


def biot_elastic_coefficients(phi: float, cfg: SchakelConfig) -> Tuple[float, float, float, float]:
    """Schakel Eq. (8)-(10): A, Q, R, and P=A+2G."""
    denom = cfg.K_f * (1.0 - phi - cfg.K_b / cfg.K_s) + phi * cfg.K_s
    A = (
        ((1.0 - phi) ** 2 * cfg.K_s * cfg.K_f - (1.0 - phi) * cfg.K_b * cfg.K_f + phi * cfg.K_s * cfg.K_b)
        / denom
        - 2.0 * cfg.G / 3.0
    )
    Q = phi * (cfg.K_s * (1.0 - phi) - cfg.K_b) * cfg.K_f / denom
    R = phi**2 * cfg.K_s * cfg.K_f / denom
    P = A + 2.0 * cfg.G
    return A, Q, R, P


def wave_slownesses(
    phi: float,
    k0: float,
    alpha_inf: float,
    cH_molL: float,
    omega: float,
    cfg: SchakelConfig,
    C_override_molL: float | None = None,
) -> Dict[str, complex]:
    """Schakel Eq. (24)-(39) slownesses and amplitude ratios."""
    dyn = dynamic_coefficients(phi, k0, alpha_inf, cH_molL, omega, cfg, C_override_molL=C_override_molL)
    k_dyn = dyn["k_dyn"]
    L = dyn["L"]
    eps_bar = dyn["eps_bar"]
    A, Q, R, P = biot_elastic_coefficients(phi, cfg)

    rho12 = phi * cfg.rho_f * (1.0 + 1j * phi * cfg.eta / (omega * cfg.rho_f * k_dyn))
    rho11 = (1.0 - phi) * cfg.rho_s - rho12
    rho22 = phi * cfg.rho_f - rho12

    E_K = cfg.eta**2 * phi**2 * L**2 / (k_dyn**2 * eps_bar * omega**2)
    rb11 = rho11 - E_K
    rb12 = rho12 + E_K
    rb22 = rho22 - E_K

    d0 = rb11 * rb22 - rb12**2
    d1 = -(P * rb22 + R * rb11 - 2.0 * Q * rb12)
    d2 = P * R - Q**2
    disc = (d1 / d2) ** 2 - 4.0 * d0 / d2
    roots_long = [(-d1 / d2 + np.sqrt(disc)) / 2.0, (-d1 / d2 - np.sqrt(disc)) / 2.0]
    roots_long = sorted(roots_long, key=lambda x: abs(x))
    s2_Pf, s2_Ps = roots_long[0], roots_long[1]

    d0_t = cfg.mu0 * eps_bar * (rb11 * rb22 - rb12**2) / cfg.G
    d1_t = -cfg.mu0 * eps_bar * rb22 - (rho11 * rho22 - rho12**2) / cfg.G
    d2_t = rho22
    disc_t = (d1_t / d2_t) ** 2 - 4.0 * d0_t / d2_t
    roots_trans = [(-d1_t / d2_t + np.sqrt(disc_t)) / 2.0, (-d1_t / d2_t - np.sqrt(disc_t)) / 2.0]
    roots_trans = sorted(roots_trans, key=lambda x: abs(x))
    s2_TM, s2_SV = roots_trans[0], roots_trans[1]

    beta_Pf = (rb11 - P * s2_Pf) / (Q * s2_Pf - rb12)
    beta_Ps = (rb11 - P * s2_Ps) / (Q * s2_Ps - rb12)
    beta_TM = (cfg.G * s2_TM - (1.0 - phi) * cfg.rho_s) / (phi * cfg.rho_f)
    beta_SV = (cfg.G * s2_SV - (1.0 - phi) * cfg.rho_s) / (phi * cfg.rho_f)

    alpha_Pf = cfg.eta * phi * L / (k_dyn * eps_bar) * (1.0 - beta_Pf)
    alpha_Ps = cfg.eta * phi * L / (k_dyn * eps_bar) * (1.0 - beta_Ps)
    alpha_TM = cfg.mu0 * cfg.eta * phi * L / (k_dyn * (cfg.mu0 * eps_bar - s2_TM)) * (1.0 - beta_TM)
    alpha_SV = cfg.mu0 * cfg.eta * phi * L / (k_dyn * (cfg.mu0 * eps_bar - s2_SV)) * (1.0 - beta_SV)

    return {
        **dyn,
        "A": A,
        "Q": Q,
        "R": R,
        "P": P,
        "rho11": rho11,
        "rho12": rho12,
        "rho22": rho22,
        "rb11": rb11,
        "rb12": rb12,
        "rb22": rb22,
        "s2_Pf": s2_Pf,
        "s2_Ps": s2_Ps,
        "s2_TM": s2_TM,
        "s2_SV": s2_SV,
        "beta_Pf": beta_Pf,
        "beta_Ps": beta_Ps,
        "beta_TM": beta_TM,
        "beta_SV": beta_SV,
        "alpha_Pf": alpha_Pf,
        "alpha_Ps": alpha_Ps,
        "alpha_TM": alpha_TM,
        "alpha_SV": alpha_SV,
    }


def se_coefficients(
    phi: float,
    k0_m2: float,
    alpha_inf: float,
    cH_molL: float,
    omega: float,
    theta_deg: float | None,
    cfg: SchakelConfig,
    kx_override: float | None = None,
    C_override_molL: float | None = None,
) -> Dict[str, complex]:
    """Solve Schakel Appendix B Eq. (B1)-(B7)."""
    state = wave_slownesses(phi, k0_m2, alpha_inf, cH_molL, omega, cfg, C_override_molL=C_override_molL)
    c_fl = math.sqrt(cfg.K_fl / cfg.rho_fl)
    if kx_override is None:
        theta = math.radians(theta_deg if theta_deg is not None else 45.0)
        k1 = omega / c_fl * math.sin(theta)
    else:
        k1 = float(kx_override)

    k3_fl = complex_sqrt_branch((omega / c_fl) ** 2 - k1**2)
    s2_E = cfg.mu0 * cfg.eps0 * cfg.eps_fl - 1j * cfg.mu0 * cfg.sigma_fl / omega
    k3_E = complex_sqrt_branch(omega**2 * s2_E - k1**2)
    k3_Pf = complex_sqrt_branch(omega**2 * state["s2_Pf"] - k1**2)
    k3_Ps = complex_sqrt_branch(omega**2 * state["s2_Ps"] - k1**2)
    k3_TM = complex_sqrt_branch(omega**2 * state["s2_TM"] - k1**2)
    k3_SV = complex_sqrt_branch(omega**2 * state["s2_SV"] - k1**2)

    Q, R, P = state["Q"], state["R"], state["P"]
    beta_Pf, beta_Ps = state["beta_Pf"], state["beta_Ps"]
    beta_TM, beta_SV = state["beta_TM"], state["beta_SV"]
    alpha_Pf, alpha_Ps = state["alpha_Pf"], state["alpha_Ps"]
    alpha_TM, alpha_SV = state["alpha_TM"], state["alpha_SV"]
    s2_Pf, s2_Ps, s2_TM, s2_SV = state["s2_Pf"], state["s2_Ps"], state["s2_TM"], state["s2_SV"]

    N1 = P - Q * (1.0 - phi) / phi + (Q - R * (1.0 - phi) / phi) * beta_Pf
    N2 = P - Q * (1.0 - phi) / phi + (Q - R * (1.0 - phi) / phi) * beta_Ps

    A_mat = np.zeros((6, 6), dtype=complex)
    A_mat[0, 1] = k3_fl
    A_mat[0, 2] = k3_Pf * (1.0 - phi + phi * beta_Pf)
    A_mat[0, 3] = k3_Ps * (1.0 - phi + phi * beta_Ps)
    A_mat[0, 4] = k1 * (1.0 - phi + phi * beta_TM)
    A_mat[0, 5] = k1 * (1.0 - phi + phi * beta_SV)

    A_mat[1, 1] = -phi * cfg.rho_fl
    A_mat[1, 2] = (Q + R * beta_Pf) * s2_Pf
    A_mat[1, 3] = (Q + R * beta_Ps) * s2_Ps

    A_mat[2, 2] = k1 * k3_Pf
    A_mat[2, 3] = k1 * k3_Ps
    A_mat[2, 4] = k1**2 - 0.5 * omega**2 * s2_TM
    A_mat[2, 5] = k1**2 - 0.5 * omega**2 * s2_SV

    A_mat[3, 2] = k1**2 - omega**2 * s2_Pf * N1 / (2.0 * cfg.G)
    A_mat[3, 3] = k1**2 - omega**2 * s2_Ps * N2 / (2.0 * cfg.G)
    A_mat[3, 4] = -k1 * k3_TM
    A_mat[3, 5] = -k1 * k3_SV

    A_mat[4, 0] = -s2_E / cfg.mu0
    A_mat[4, 4] = alpha_TM * s2_TM / cfg.mu0
    A_mat[4, 5] = alpha_SV * s2_SV / cfg.mu0

    A_mat[5, 0] = -k3_E
    A_mat[5, 2] = k1 * alpha_Pf
    A_mat[5, 3] = k1 * alpha_Ps
    A_mat[5, 4] = -k3_TM * alpha_TM
    A_mat[5, 5] = -k3_SV * alpha_SV

    b_vec = np.array([k3_fl, phi * cfg.rho_fl, 0.0, 0.0, 0.0, 0.0], dtype=complex)
    try:
        x = np.linalg.solve(A_mat, b_vec)
    except np.linalg.LinAlgError:
        x = np.full(6, np.nan + 1j * np.nan, dtype=complex)

    return {
        **state,
        "k1": k1,
        "k3_fl": k3_fl,
        "k3_E": k3_E,
        "k3_Pf": k3_Pf,
        "k3_Ps": k3_Ps,
        "k3_TM": k3_TM,
        "k3_SV": k3_SV,
        "s2_E": s2_E,
        "R_E": x[0],
        "R_M": x[1],
        "T_Pf": x[2],
        "T_Ps": x[3],
        "T_TM": x[4],
        "T_SV": x[5],
        "matrix_cond": np.linalg.cond(A_mat) if np.all(np.isfinite(A_mat)) else np.nan,
    }


def incident_acoustic_flux(coeff: Mapping[str, complex], omega: float, cfg: SchakelConfig) -> float:
    """Incident fluid P-wave vertical flux for unit incident scalar potential."""
    return 0.5 * cfg.rho_fl * omega**2 * float(np.real(coeff["k3_fl"]))


def reflected_em_flux_coefficient(coeff: Mapping[str, complex], omega: float, cfg: SchakelConfig) -> float:
    pin = incident_acoustic_flux(coeff, omega, cfg)
    r_e = coeff["R_E"]
    e1 = -coeff["k3_E"] * r_e
    h2 = -(coeff["s2_E"] / cfg.mu0) * r_e
    flux = 0.5 * np.real(e1 * np.conj(h2))
    return float(flux / pin)


def reflected_acoustic_flux_coefficient(coeff: Mapping[str, complex]) -> float:
    return float(-abs(coeff["R_M"]) ** 2)


def transmitted_tm_flux_coefficient(coeff: Mapping[str, complex], omega: float, cfg: SchakelConfig) -> float:
    """Schakel orthodox transmitted TM flux coefficient T_E^{TM,TM}."""
    pin = incident_acoustic_flux(coeff, omega, cfg)
    t_tm = coeff["T_TM"]
    k1 = coeff["k1"]
    k3_tm = coeff["k3_TM"]
    s2_tm = coeff["s2_TM"]
    alpha_tm = coeff["alpha_TM"]

    e1 = -k3_tm * alpha_tm * t_tm
    h2 = alpha_tm * s2_tm / cfg.mu0 * t_tm
    u1 = k3_tm * t_tm
    u3 = k1 * t_tm
    tau31 = -2.0 * cfg.G * (k1**2 - 0.5 * omega**2 * s2_tm) * t_tm
    tau33 = -2.0 * cfg.G * k1 * k3_tm * t_tm
    flux = 0.5 * np.real(
        e1 * np.conj(h2)
        + tau31 * np.conj(1j * omega * u1)
        + tau33 * np.conj(1j * omega * u3)
    )
    return float(flux / pin)


def strict_energy_flux_coefficients(
    phi: float,
    k0_m2: float,
    alpha_inf: float,
    cH_molL: float,
    omega: float,
    theta_deg: float,
    cfg: SchakelConfig,
    C_override_molL: float | None = None,
) -> Dict[str, float | complex]:
    coeff = se_coefficients(phi, k0_m2, alpha_inf, cH_molL, omega, theta_deg, cfg, C_override_molL=C_override_molL)
    re_ee = reflected_em_flux_coefficient(coeff, omega, cfg)
    re_pr_pr = reflected_acoustic_flux_coefficient(coeff)
    te_tm_tm = transmitted_tm_flux_coefficient(coeff, omega, cfg)
    return {
        **coeff,
        "RE_EE": re_ee,
        "TE_TM_TM": te_tm_tm,
        "minus_TE_TM_TM": -te_tm_tm,
        "RE_Pr_Pr": re_pr_pr,
        "Pin_flux": incident_acoustic_flux(coeff, omega, cfg),
    }


def reference_energy_flux_coefficients(omega: float = 1.0e6, theta_deg: float = 45.0) -> Dict[str, float | complex]:
    cfg, baseline = schakel_reference_config()
    return strict_energy_flux_coefficients(
        baseline["phi"],
        baseline["k0_m2"],
        baseline["alpha_inf"],
        baseline["cH_molL"],
        omega,
        theta_deg,
        cfg,
        C_override_molL=baseline["C_molL"],
    )


def table_ii_validation_rows() -> list[Dict[str, float]]:
    rows = []
    computed_keys = {"RE_EE", "TE_TM_TM", "RE_Pr_Pr"}
    te_keys = [key for key in next(iter(TABLE_II_REFERENCE.values())) if key.startswith("TE_")]
    for frequency_hz, reference in TABLE_II_REFERENCE.items():
        result = reference_energy_flux_coefficients(omega=2.0 * math.pi * frequency_hz, theta_deg=30.0)
        row = {"frequency_hz": frequency_hz, "omega_rad_s": 2.0 * math.pi * frequency_hz}
        for key in reference:
            row[key] = float(result[key]) if key in computed_keys else reference[key]
            row[f"{key}_paper"] = reference[key]
            if key in computed_keys:
                row[f"{key}_computed"] = float(result[key])
                row[f"{key}_relative_error"] = (float(result[key]) - reference[key]) / reference[key]
                row[f"{key}_source"] = "computed"
            else:
                row[f"{key}_computed"] = np.nan
                row[f"{key}_relative_error"] = np.nan
                row[f"{key}_source"] = "paper_table"
        row["energy_balance"] = float(sum(row[key] for key in te_keys) - row["RE_EE"] - row["RE_Pr_Pr"])
        row["energy_balance_residual"] = row["energy_balance"] - 1.0
        row["note"] = "Core Schakel formula checks: RE_EE, TE_TM_TM, RE_Pr_Pr. Other TE entries reproduce Table II values."
        rows.append(row)
    return rows


def table_iii_validation_rows() -> list[Dict[str, float]]:
    rows = []
    computed_keys = {"RE_EE", "TE_TM_TM", "RE_Pr_Pr"}
    te_keys = [key for key in next(iter(TABLE_III_REFERENCE.values())) if key.startswith("TE_")]
    for theta_deg, reference in TABLE_III_REFERENCE.items():
        result = reference_energy_flux_coefficients(omega=1.0e6, theta_deg=theta_deg)
        row = {"omega_rad_s": 1.0e6, "theta_deg": theta_deg}
        for key in reference:
            row[key] = float(result[key]) if key in computed_keys else reference[key]
            row[f"{key}_paper"] = reference[key]
            if key in computed_keys:
                row[f"{key}_computed"] = float(result[key])
                row[f"{key}_relative_error"] = (float(result[key]) - reference[key]) / reference[key]
                row[f"{key}_source"] = "computed"
            else:
                row[f"{key}_computed"] = np.nan
                row[f"{key}_relative_error"] = np.nan
                row[f"{key}_source"] = "paper_table"
        row["energy_balance"] = float(sum(row[key] for key in te_keys) - row["RE_EE"] - row["RE_Pr_Pr"])
        row["energy_balance_residual"] = row["energy_balance"] - 1.0
        row["note"] = "Core Schakel formula checks: RE_EE, TE_TM_TM, RE_Pr_Pr. Other TE entries reproduce Table III values."
        rows.append(row)
    return rows


def figure2_angle_scan(n: int = 201) -> pd.DataFrame:
    cfg, baseline = schakel_reference_config()
    omega = 2.0 * math.pi * 10.0
    c_fl = math.sqrt(cfg.K_fl / cfg.rho_fl)
    rows = []
    t_tm_values = []
    for sin_theta in np.linspace(0.0, 2.0, int(n)):
        coeff = se_coefficients(
            baseline["phi"],
            baseline["k0_m2"],
            baseline["alpha_inf"],
            baseline["cH_molL"],
            omega,
            None,
            cfg,
            kx_override=omega / c_fl * float(sin_theta),
            C_override_molL=baseline["C_molL"],
        )
        t_tm_values.append(coeff["T_TM"])
        rows.append(
            {
                "sin_theta": float(sin_theta),
                "frequency_hz": 10.0,
                "R_E_real": float(np.real(coeff["R_E"])),
                "R_E_imag": float(np.imag(coeff["R_E"])),
                "R_E_abs": float(abs(coeff["R_E"])),
                "R_E_phase_rad": float(np.angle(coeff["R_E"])),
                "T_TM_real": float(np.real(coeff["T_TM"])),
                "T_TM_imag": float(np.imag(coeff["T_TM"])),
                "T_TM_abs": float(abs(coeff["T_TM"])),
                "T_TM_phase_rad": float(np.angle(coeff["T_TM"])),
                "matrix_cond": float(coeff["matrix_cond"]),
            }
        )
    df = pd.DataFrame(rows)
    t_phase = np.unwrap(np.angle(np.asarray(t_tm_values, dtype=complex)))
    t_phase = np.where(t_phase > 4.0 * math.pi, t_phase - 2.0 * math.pi, t_phase)
    t_phase = np.where(t_phase < 0.0, t_phase + 2.0 * math.pi, t_phase)
    df["T_TM_phase_rad"] = t_phase
    return df


def plot_figure2(df: pd.DataFrame, outpath: Path) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(7.6, 5.6), constrained_layout=True)
    axes[0, 0].plot(df["sin_theta"], df["R_E_abs"], color="black", linewidth=1.0)
    axes[0, 1].plot(df["sin_theta"], df["R_E_phase_rad"], color="black", linewidth=1.0)
    axes[1, 0].plot(df["sin_theta"], df["T_TM_abs"], color="black", linewidth=1.0)
    axes[1, 1].plot(df["sin_theta"], df["T_TM_phase_rad"], color="black", linewidth=1.0)
    axes[0, 0].set_ylabel(r"$|R^E|$ (V/m$^2$)")
    axes[0, 1].set_ylabel(r"$\angle R^E$")
    axes[1, 0].set_ylabel(r"$|T^{TM}|$")
    axes[1, 1].set_ylabel(r"$\angle T^{TM}$")
    for ax in axes.ravel():
        ax.set_xlabel(r"$\sin\theta$")
        ax.set_xlim(0.0, 2.0)
        ax.grid(False)
    axes[0, 1].set_ylim(-math.pi, math.pi)
    axes[0, 0].set_ylim(0.0, 20.0)
    axes[1, 0].set_ylim(0.0, 1.0e-9)
    axes[1, 1].set_ylim(0.0, 4.0 * math.pi)
    fig.suptitle("Schakel & Smeulders Fig. 2 reproduction")
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def scan_definitions() -> Dict[str, np.ndarray]:
    return {
        "electrolyte_concentration_molL": np.logspace(-6, 0, 100),
        "viscosity_Pa_s": np.logspace(-6, 0, 100),
        "permeability_m2": np.logspace(-16, -10, 120),
        "tortuosity": np.linspace(1.0, 5.0, 100),
        "pH": np.linspace(2.0, 10.0, 100),
        "pore_fluid_bulk_modulus_Pa": np.logspace(7, 11, 100),
        "porosity": np.linspace(0.10, 0.85, 100),
    }


def run_single_scan(scan_name: str, values: Sequence[float]) -> pd.DataFrame:
    cfg, baseline = schakel_reference_config()
    rows = []
    for value in values:
        local_cfg = replace(cfg)
        phi = baseline["phi"]
        k0_m2 = baseline["k0_m2"]
        alpha_inf = baseline["alpha_inf"]
        cH = baseline["cH_molL"]
        c_override = baseline["C_molL"]
        if scan_name == "electrolyte_concentration_molL":
            c_override = float(value)
        elif scan_name == "viscosity_Pa_s":
            local_cfg.eta = float(value)
        elif scan_name == "permeability_m2":
            k0_m2 = float(value)
        elif scan_name == "tortuosity":
            alpha_inf = float(value)
        elif scan_name == "pH":
            cH = h_concentration_for_ph(float(value))
        elif scan_name == "pore_fluid_bulk_modulus_Pa":
            local_cfg.K_f = float(value)
        elif scan_name == "porosity":
            phi = float(value)
        else:
            raise ValueError(f"Unknown scan: {scan_name}")

        coeff = strict_energy_flux_coefficients(
            phi,
            k0_m2,
            alpha_inf,
            cH,
            baseline["omega"],
            baseline["theta_deg"],
            local_cfg,
            C_override_molL=c_override,
        )
        rows.append(
            {
                "scan": scan_name,
                "scan_value": float(value),
                "Porosity_used": phi,
                "k0_m2": k0_m2,
                "Permeability_mD": k0_m2 / 9.869233e-16,
                "Tortuosity": alpha_inf,
                "C_molL": float(coeff["C_molL"]),
                "pH": float(coeff["pH"]),
                "eta_Pa_s": local_cfg.eta,
                "K_f_Pa": local_cfg.K_f,
                "Lambda_m": float(coeff["Lambda"]),
                "debye_d_m": float(coeff["debye_d"]),
                "L_abs": abs(coeff["L"]),
                "R_E_abs": abs(coeff["R_E"]),
                "T_TM_abs": abs(coeff["T_TM"]),
                "RE_EE": coeff["RE_EE"],
                "TE_TM_TM": coeff["TE_TM_TM"],
                "minus_TE_TM_TM": coeff["minus_TE_TM_TM"],
                "RE_Pr_Pr": coeff["RE_Pr_Pr"],
                "matrix_cond": float(coeff["matrix_cond"]),
            }
        )
    return pd.DataFrame(rows)


def plot_two_by_two(
    df: pd.DataFrame,
    scan_names: Sequence[str],
    outpath: Path,
    title: str,
    xlabels: Mapping[str, str],
    xscale: Mapping[str, str],
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(8.6, 6.3), constrained_layout=True)
    metrics = [("RE_EE", r"$R_E^{E,E}$"), ("minus_TE_TM_TM", r"$-T_E^{TM,TM}$")]
    for col, scan_name in enumerate(scan_names):
        part = df[df["scan"] == scan_name]
        for row, (metric, ylabel) in enumerate(metrics):
            ax = axes[row, col]
            ax.plot(part["scan_value"], part[metric], color="black", linewidth=1.3)
            ax.set_xlabel(xlabels[scan_name])
            ax.set_ylabel(ylabel)
            if xscale.get(scan_name) == "log":
                ax.set_xscale("log")
            ax.grid(False)
    fig.suptitle(title)
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def plot_porosity(df: pd.DataFrame, outpath: Path) -> None:
    part = df[df["scan"] == "porosity"]
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.5), constrained_layout=True)
    axes[0].plot(part["scan_value"], part["RE_EE"], color="black", linewidth=1.3)
    axes[1].plot(part["scan_value"], part["minus_TE_TM_TM"], color="black", linewidth=1.3)
    axes[0].set_ylabel(r"$R_E^{E,E}$")
    axes[1].set_ylabel(r"$-T_E^{TM,TM}$")
    for ax in axes:
        ax.set_xlabel("porosity")
    fig.suptitle("Schakel Fig. 7 reproduction")
    fig.savefig(outpath, dpi=300)
    plt.close(fig)


def run_strict_sensitivity(outdir: Path) -> pd.DataFrame:
    outdir.mkdir(parents=True, exist_ok=True)
    fig2 = figure2_angle_scan()
    fig2.to_csv(outdir / "schakel2010_fig2_angle_scan.csv", index=False)
    plot_figure2(fig2, outdir / "schakel2010_fig2_reproduction.png")

    pd.DataFrame(table_ii_validation_rows()).to_csv(outdir / "schakel2010_tableII_validation.csv", index=False)
    pd.DataFrame(table_iii_validation_rows()).to_csv(outdir / "schakel2010_tableIII_validation.csv", index=False)

    frames = [run_single_scan(name, values) for name, values in scan_definitions().items()]
    df = pd.concat(frames, ignore_index=True)
    df.to_csv(outdir / "schakel2010_strict_flux_coefficients.csv", index=False)

    summary_rows = []
    for scan, part in df.groupby("scan"):
        for metric in ["RE_EE", "minus_TE_TM_TM", "L_abs"]:
            idx = part[metric].idxmax()
            min_val = part[metric].min()
            max_val = part[metric].max()
            summary_rows.append(
                {
                    "scan": scan,
                    "metric": metric,
                    "min": min_val,
                    "max": max_val,
                    "max_over_min": max_val / min_val if min_val > 0 else np.nan,
                    "scan_value_at_max": df.loc[idx, "scan_value"],
                }
            )
    pd.DataFrame(summary_rows).to_csv(outdir / "schakel2010_strict_flux_summary.csv", index=False)

    xlabels = {
        "electrolyte_concentration_molL": "C (mol/L)",
        "viscosity_Pa_s": r"$\eta$ (Pa s)",
        "permeability_m2": r"$k_0$ (m$^2$)",
        "tortuosity": r"$\alpha_\infty$",
        "pH": "pH",
        "pore_fluid_bulk_modulus_Pa": r"$K_f$ (Pa)",
    }
    xscale = {
        "electrolyte_concentration_molL": "log",
        "viscosity_Pa_s": "log",
        "permeability_m2": "log",
        "pore_fluid_bulk_modulus_Pa": "log",
    }
    plot_two_by_two(
        df,
        ["electrolyte_concentration_molL", "viscosity_Pa_s"],
        outdir / "fig4_strict_concentration_viscosity.png",
        "Schakel Fig. 4 reproduction",
        xlabels,
        xscale,
    )
    plot_two_by_two(
        df,
        ["permeability_m2", "tortuosity"],
        outdir / "fig5_strict_permeability_tortuosity.png",
        "Schakel Fig. 5 reproduction",
        xlabels,
        xscale,
    )
    plot_two_by_two(
        df,
        ["pH", "pore_fluid_bulk_modulus_Pa"],
        outdir / "fig6_strict_ph_bulk_modulus.png",
        "Schakel Fig. 6 reproduction",
        xlabels,
        xscale,
    )
    plot_porosity(df, outdir / "fig7_strict_porosity.png")
    return df


def write_readme(outdir: Path) -> None:
    text = """# Strict Schakel & Smeulders (2010) Reproduction

This directory is generated by `schakel2010_strict_sensitivity.py`.

The script is a standalone Schakel & Smeulders (2010) reproduction. It uses
the paper's Table I parameters, Appendix A dynamic coefficients, Appendix B
boundary-value system, Fig. 2 potential-coefficient definitions, and Eq.
(47)-(51) vertical energy-flux definitions.

For Table II and Table III, the independently recomputed columns are
`RE_EE_computed`, `TE_TM_TM_computed`, and `RE_Pr_Pr_computed`; these are
compared with the paper values and relative errors are reported. The remaining
orthodox/interference transmission entries are included as the paper table
values so the complete published tables and Eq. (51) balances are present
without presenting an unverified cross-flux derivation as computed output.

Outputs:

- `schakel2010_fig2_angle_scan.csv` and `schakel2010_fig2_reproduction.png`.
- `schakel2010_tableII_validation.csv`.
- `schakel2010_tableIII_validation.csv`.
- `fig4_strict_concentration_viscosity.png`.
- `fig5_strict_permeability_tortuosity.png`.
- `fig6_strict_ph_bulk_modulus.png`.
- `fig7_strict_porosity.png`.
- `schakel2010_strict_flux_coefficients.csv`.
- `schakel2010_strict_flux_summary.csv`.
"""
    (outdir / "README.md").write_text(text, encoding="utf-8")

    cn = """# 严格 Schakel & Smeulders (2010) 复现说明

本目录由 `schakel2010_strict_sensitivity.py` 生成，只作为 Schakel &
Smeulders (2010) 原文公式和参数的独立复现。

实现范围包括：

- Table I 参数。
- Appendix A 的动态渗透率、动电耦合系数和动态电导率。
- Eq. (24)-(39) 的慢度与幅值比。
- Appendix B 的六阶界面边界值问题。
- Fig. 2 中 `R_E` 与 `T_TM` 势函数系数随 `sin(theta)` 的幅值和相位。
- Eq. (47)-(51) 的垂向能流系数，并输出 Table II 与 Table III 数值验证。

其中 `R_E` 和 `T_TM` 是 Schakel Eq. (45) 定义的反射 EM 矢量势系数和透射
TM 矢量势系数；能流表使用 Schakel Eq. (47)-(51) 的 orthodox/interference
flux 定义。

Table II 和 Table III 中，本脚本独立重算并校验 `RE_EE_computed`、
`TE_TM_TM_computed` 和 `RE_Pr_Pr_computed`；其余 orthodox/interference
transmission 项作为原文表格数值一起输出，用于完整保留论文表格和 Eq. (51)
能量守恒检查，而不把未完成复核的 cross-flux 推导冒充为已计算结果。
"""
    (outdir / "结果说明.md").write_text(cn, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", default="sensitivity_results/schakel2010_strict_reproduction")
    args = parser.parse_args()
    outdir = Path(args.outdir)
    run_strict_sensitivity(outdir)
    write_readme(outdir)


if __name__ == "__main__":
    main()
