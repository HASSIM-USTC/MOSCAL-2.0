import time
import numpy as np
import subprocess
import itertools
import json
from scipy.linalg import expm as matrix_exp
import sys
from deom import benchmark, convert, complex_2_json, init_qmd, init_qmd_quad


for npsd, lmax, theta, eta, init, ferr, beta in itertools.product(
    [1],  # npsd
    [25],  # lmax
    [4 / 3],  # theta
    [0.25],  # eta
    [1],  # init, 1 for absorption, 0 for emission
    [1e-16],  # ferr
    [1]  # beta
):
    nmod = 1
    weg = 1
    omgb = 1
    ti_gam = 15
    ti_eta = 10
    sdip_val, bdip0_val, bdip1_val, bdip2_val = 1, 1, 1, 1
    dt = 0.001

    if (init == 1):
        omgb_bo = 1
    if (init == 0):
        omgb_bo = theta

    if (init == 1):
        alp0 = eta * theta * theta
        alp1 = -np.sqrt(2 * eta * omgb) * theta * theta
        alp2 = omgb * (theta * theta - 1) / 2
    if (init == 0):
        alp0 = eta
        alp1 = np.sqrt(2 * eta * omgb / theta)
        alp2 = -omgb * (theta * theta - 1) / 2 / theta
    expn = np.array([
        0.3477781831017336 - 0.963036125471659j,
        0.3477781831016224 + 0.9630361254716366j
    ], dtype=complex)
    etal = np.array([
        0.2703547700830425 - 0.167441733506598j,
        0.8145601320944829 + 0.167441733506598j
    ], dtype=complex)
    etar = np.array([
        0.8145601320944829 - 0.167441733506598j,
        0.2703547700830425 + 0.167441733506598j
    ], dtype=complex)
    etaa = np.array([0.3180069744932481, 0.831591812680642])

    mode = np.zeros_like(expn, dtype=int)

    hams = np.zeros((2, 2), dtype=complex)
    hams[0, 0] = 1

    mat_equ = np.zeros((2, 2), dtype=complex)
    mat_equ = matrix_exp(- beta * hams) / np.trace(matrix_exp(- beta * hams))

    qmds = np.zeros((nmod, 2, 2), dtype=complex)
    if (init == 1):
        qmds[0, 0, 0] = 1
    if (init == 0):
        qmds[0, 1, 1] = 1

    qmd2 = np.zeros((nmod, nmod, 2, 2), dtype=complex)
    for i in range(nmod):
        for j in range(nmod):
            qmd2[i, j, :, :] = qmds[i, :, :]

    trun = 0
    nmax = 2000000

    sdip = np.zeros((nmod, 2, 2), dtype=float)
    sdip[0, 0, 1] = sdip_val * 1.0
    sdip[0, 1, 0] = sdip_val * 1.0

    pdip0 = np.zeros((nmod, 2, 2), dtype=float)
    pdip1 = np.zeros((nmod, 2, 2), dtype=float)
    pdip2 = np.zeros((nmod, 2, 2), dtype=float)

    bdip0 = bdip0_val * np.zeros(12 + 2 * npsd, dtype=float)
    bdip1 = bdip1_val * np.zeros(12 + 2 * npsd, dtype=float)
    bdip2 = bdip2_val * np.zeros(12 + 2 * npsd, dtype=float)

    sdipn = np.zeros((nmod, 2, 2), dtype=float)
    pdipn = np.zeros((nmod, 2, 2), dtype=float)
    bdipn = np.zeros(12 + 2 * npsd, dtype=float)

    # proprho
    json_init = {
        "nmax": nmax,
        "lmax": lmax,
        "alp0": alp0,
        "alp1": alp1,
        "alp2": alp2,
        "ferr": ferr,
        "nind": len(expn),
        "nmod": nmod,
        "inistate": 1,
        "imag": True,
        "equilibrium": {
            "ti": 0,
            "tf": beta,
            "dt": dt,
        },
        "expn": complex_2_json(expn),
        "read_rho0": True,
        "rho0": complex_2_json(mat_equ),
        "ham1": complex_2_json(hams),
        "coef_abs": complex_2_json(etaa),
        "renormalize": complex_2_json((alp0 + alp2 * np.sum(etal)) * qmds),
    }

    init_qmd(json_init, qmds, qmds, mode, 2, etaa, etal, etar)
    init_qmd_quad(json_init, qmd2, qmd2, qmd2, mode,
                  2, len(expn), nmod, etaa, etal, etar)

    magic_str = '{}-imag'.format(ferr)
    with open('input.json', 'w') as f:
        json.dump(json_init, f, indent=4, default=convert)
    cmd = r'export OMP_NUM_THREADS={}'.format(
        8 if len(sys.argv) == 1 else sys.argv[1])
    cmd += '&&' + r'/usr/local/bin/jemalloc.sh  ../thermo_bose_quad_2.out'
    start_time = time.time()
    with open('out-{}'.format(magic_str), "w") as outfile:
        result = subprocess.call(cmd, shell=True, stdout=outfile)
    np.savetxt('time-{}'.format(magic_str), [time.time() - start_time])
    benchmark('prop-rho-eq1.dat', magic_str, 'bose_quad')
