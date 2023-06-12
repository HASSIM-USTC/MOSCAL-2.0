import numpy as np
import sympy as sp
import itertools
import json
from deom import benchmark, convert, decompose_spe, complex_2_json
import sys
import time
import subprocess


def mea(omgm, gamm, temp):
    coef = np.zeros(8, dtype=complex)
    expn = np.zeros(2, dtype=complex)
    etal = np.zeros(2, dtype=complex)
    etar = np.zeros(2, dtype=complex)
    etaa = np.zeros(2, dtype=float)

    # para detail in Minimum-exponents ansatz for molecular dynamics and quantum dissipation
    if gamm * gamm - 4 * omgm * omgm > 0:
        gamp = 1. / 2. * (gamm + np.sqrt(gamm * gamm - 4 * omgm * omgm + 0.j))
        gamn = 1. / 2. * (gamm - np.sqrt(gamm * gamm - 4 * omgm * omgm + 0.j))
        expn[0] = gamp
        expn[1] = gamn
        etal[0] = -omgm / (gamp - gamn) * (1 * temp / gamp - 1.j / 2.)
        etal[1] = omgm / (gamp - gamn) * (1 * temp / gamn - 1.j / 2.)
        etar[0] = etal[0].conj()
        etar[1] = etal[1].conj()
        etaa[0] = abs(etal[0])
        etaa[1] = abs(etal[1])
        coef[1] = -1.j * omgm / 2. * (2. * etal[0] - 2. * etar[0] -
                                      2. * etal[1] + 2. * etar[1])
        coef[2] = -1.j * omgm / 2. * (
            -2. * etal[1] * gamn * gamn / omgm / omgm + 2. * etar[1] * gamn *
            gamn / omgm / omgm + 2. * etal[0] - 2. * etar[0])
        coef[3] = -1.j * omgm / 2. * (
            -2. * etal[0] * gamp * gamp / omgm / omgm + 2. * etar[0] * gamp *
            gamp / omgm / omgm + 2. * etal[1] - 2. * etar[1])
    else:
        gamp = 1. / 2. * (gamm + np.sqrt(gamm * gamm - 4 * omgm * omgm + 0.j))
        gamn = 1. / 2. * (gamm - np.sqrt(gamm * gamm - 4 * omgm * omgm + 0.j))
        expn[0] = gamp
        expn[1] = gamn
        etal[0] = -omgm / (gamp - gamn) * (1 * temp / gamp - 1.j / 2.)
        etal[1] = omgm / (gamp - gamn) * (1 * temp / gamn - 1.j / 2.)
        etar[1] = etal[0].conj()
        etar[0] = etal[1].conj()
        etaa[0] = abs(etal[0])
        etaa[1] = abs(etal[1])
        coef[1] = -1.j * omgm / 2. * (2. * etal[0] - 2. * etar[0] -
                                      2. * etal[1] + 2. * etar[1])
        coef[2] = -1.j * omgm / 2. * (
            -2. * etal[1] * gamn * gamn / omgm / omgm + 2. * etar[1] * gamn *
            gamn / omgm / omgm + 2. * etal[0] - 2. * etar[0])
        coef[3] = -1.j * omgm / 2. * (
            -2. * etal[0] * gamp * gamp / omgm / omgm + 2. * etar[0] * gamp *
            gamp / omgm / omgm + 2. * etal[1] - 2. * etar[1])

    # return
    return etal, etar, etaa, expn, coef


for npsd, lmax, (lmax_fp, lmax_ma), (theta, eta), (
        sdip_val, bdip0_val, bdip1_val, bdip2_val
), init, dt in itertools.product(
    [1],  # npsd
    [100],  # lmax
    [(3, 8)],
    [(4 / 3, 0.25)],  # theta, eta
    [(1, 1, 1, 1)],  # bdip1_val, bdip2_val
    [1],  # init, 1 for absorption, 0 for emission
    [0.01]  # dt
):
    nmod = 1

    weg = 1
    temp = 1
    beta = 1 / temp
    omgb = 1
    ti_gam = 15
    ti_eta = 10

    if (init == 1):
        omgb_BO = 1
    if (init == 0):
        omgb_BO = theta

    w_sp, lamd_sp, gamd_sp, beta_sp = sp.symbols(
        r"\omega , \lambda, \gamma, \beta", real=True)

    if (init == 1):
        alp0 = eta * theta * theta
        alp1 = -np.sqrt(2 * eta * omgb) * theta * theta
        alp2 = omgb * (theta * theta - 1) / 2
    if (init == 0):
        alp0 = eta
        alp1 = np.sqrt(2 * eta * omgb / theta)
        alp2 = -omgb * (theta * theta - 1) / 2 / theta

    phixx_sp = lamd_sp * gamd_sp * w_sp / (w_sp**2 + gamd_sp**2)
    spe_vib_sp = phixx_sp * omgb / omgb_BO
    sp_para_dict = {lamd_sp: ti_eta, gamd_sp: ti_gam}

    condition_dict = {}
    para_dict = {'beta': beta}

    etal1, etar1, etaa1, expn1, coef = mea(
        omgb_BO, omgb * ti_eta / ti_gam, temp)
    etal2, etar2, etaa2, expn2 = decompose_spe(
        spe_vib_sp, w_sp, sp_para_dict, para_dict, condition_dict, npsd)
    coef[0] = omgb_BO
    etal = np.append(etal1, etal2)
    etar = np.append(etar1, etar2)
    etaa = np.append(np.ones_like(etaa1), etaa2)
    expn = np.append(np.zeros_like(expn1), expn2)
    mode = np.zeros_like(expn, dtype=int)
    beta = para_dict['beta']

    hams = np.zeros((2, 2), dtype=complex)
    hams[0, 0] = 1

    qmds = np.zeros((nmod, 2, 2), dtype=complex)
    if (init == 1):
        qmds[0, 0, 0] = 1
    if (init == 0):
        qmds[0, 1, 1] = 1

    trun = 0
    nmax = 1000000
    ferr = 1e-6

    sdip = np.zeros((nmod, 2, 2), dtype=float)
    sdip[0, 0, 1] = sdip_val * 1.0
    sdip[0, 1, 0] = sdip_val * 1.0

    pdip0 = np.zeros((nmod, 2, 2), dtype=float)
    pdip1 = np.zeros((nmod, 2, 2), dtype=float)
    pdip2 = np.zeros((nmod, 2, 2), dtype=float)
    pdip0[0, 0, 0] = 1.0
    pdip0[0, 1, 1] = 1.0
    pdip1[0, 0, 0] = 1.0
    pdip1[0, 1, 1] = 1.0
    pdip2[0, 0, 0] = 1.0
    pdip2[0, 1, 1] = 1.0

    bdip0 = bdip0_val * np.ones(12 + 2 * npsd, dtype=float)
    bdip1 = bdip1_val * np.ones(12 + 2 * npsd, dtype=float)
    bdip2 = bdip2_val * np.ones(12 + 2 * npsd, dtype=float)

    json_init = {
        "nmax": nmax,
        "lmax": lmax,
        "lmax_fp": lmax_fp,
        "lmax_ma": lmax_ma,
        "ferr": ferr,
        "filter": True,
        "nind": len(expn),
        "nmod": nmod,
        "inistate": init,
        "alp0": alp0,
        "alp1": alp1,
        "alp2": alp2,
        "equilibrium": {
            "sc2": False,
            "dt-method": True,
            "ti": 0,
            "tf": 10,
            "dt": dt,
            "backup": True,
        },
        "coef0": complex_2_json(coef[1]),
        "coef1": complex_2_json(coef[2]),
        "coef2": complex_2_json(coef[3]),
        "coef3": complex_2_json(ti_eta / 2.0 * omgb / omgb_BO),
        "expn": complex_2_json(expn),
        "ham1": complex_2_json(hams),
        "qmd1": complex_2_json(qmds),
        "coef_lft": complex_2_json(etal),
        "coef_rht": complex_2_json(etar),
        "coef_abs": complex_2_json(etaa),
        "spectrum": True,
        "spectrum-data": {
            "dipole": {
                "sdip_cub": complex_2_json(sdip),
                "bdip0_cub": complex_2_json(bdip0),
                "bdip1_cub": complex_2_json(bdip1),
                "bdip2_cub": complex_2_json(bdip2),
                "pdip0_cub": complex_2_json(pdip0),
                "pdip1_cub": complex_2_json(pdip1),
                "pdip2_cub": complex_2_json(pdip2),
            },
            "dipole1": {
                "sdip_cub": complex_2_json(sdip),
                "bdip0_cub": complex_2_json(bdip0),
                "bdip1_cub": complex_2_json(bdip1),
                "bdip2_cub": complex_2_json(bdip2),
                "pdip0_cub": complex_2_json(pdip0),
                "pdip1_cub": complex_2_json(pdip1),
                "pdip2_cub": complex_2_json(pdip2),
            },
            "if-time": True,
            "time": {
                "ti": 0,
                "tf": 10,
                "dt": dt,
                "lcr": 'l',
                "filter_ferr": ferr,
            },
            "file": "prop-pol-1.dat",
        },
    }

    magic_str = 'corr'
    with open('input.json', 'w') as f:
        json.dump(json_init, f, indent=4, default=convert)
    cmd = r'export OMP_NUM_THREADS={}'.format(
        8 if len(sys.argv) == 1 else sys.argv[1])
    cmd += '&&' + r'/usr/local/bin/jemalloc.sh  ../bsm_2.out'
    start_time = time.time()
    with open('out-{}'.format(magic_str), "w") as outfile:
        result = subprocess.call(cmd, shell=True, stdout=outfile)
    np.savetxt('time-{}'.format(magic_str), [time.time() - start_time])
    benchmark("prop-pol-1.dat", magic_str, 'bsm')
