import numpy as np
import sympy as sp
import itertools
import json
from deom import benchmark, convert, decompose_spe, complex_2_json, init_qmd
import sys
import time
import subprocess


for npsd, lmax, ferr, (gam, eta), dt in itertools.product(
    [2],  # npsd
    [100],  # lmax
    [1e-10],  # ferr
    [(1, 1)],  # gam,eta
    [0.01]  # dt
):
    nmod = 1
    temp = 1
    beta = 1 / temp
    OMG = 50

    w_sp, eta_sp, gamma_sp, beta_sp = sp.symbols(
        r"\omega, \eta, \gamma, \beta", real=True)

    phixx_sp = 2 * eta_sp * gamma_sp / (gamma_sp - sp.I * w_sp)
    spe_vib_sp = phixx_sp
    sp_para_dict = {eta_sp: eta, gamma_sp: gam}
    condition_dict = {}
    para_dict = {'beta': beta}
    etal, etar, etaa, expn = decompose_spe(spe_vib_sp, w_sp, sp_para_dict, para_dict,
                                           condition_dict, npsd)

    mode = np.zeros_like(expn, dtype=int)
    beta = para_dict['beta']

    rho0 = np.zeros((2, 2), dtype=complex)
    rho0[0, 0] = 1
    rho0[1, 1] = 0

    hams = np.zeros((2, 2), dtype=complex)
    hams[0, 0] = 1
    hams[0, 1] = 1
    hams[1, 0] = 1
    hams[1, 1] = -1

    qmds = np.zeros((nmod, 2, 2), dtype=complex)
    qmds[0, 0, 1] = 1
    qmds[0, 1, 0] = 1

    nmax = 1000000

    twsg = np.zeros((len(expn), nmod))
    for i in range(len(expn)):
        for j in range(nmod):
            twsg[i, j] = i + j * len(expn)

    json_init = {
        "nmax": nmax,
        "lmax": lmax,
        "ferr": ferr,
        "filter": True,
        "nind": len(expn),
        "nmod": nmod,
        "equilibrium": {
            "sc2": False,
            "dt-method": True,
            "ti": 0,
            "tf": 25,
            "dt": dt,
            "backup": True,
        },
        "expn": complex_2_json(expn),
        "ham1": complex_2_json(hams),
        "coef_abs": complex_2_json(etaa),
        "read_rho0": True,
        "rho0": complex_2_json(rho0),
    }

    init_qmd(json_init, qmds, qmds, mode, 2, etaa, etal, etar)

    magic_str = '{}-{}-dt'.format(ferr, eta)
    with open('input.json', 'w') as f:
        json.dump(json_init, f, indent=4, default=convert)
    cmd = r'export OMP_NUM_THREADS={}'.format(
        8 if len(sys.argv) == 1 else sys.argv[1])
    cmd += '&&' + r'/usr/local/bin/jemalloc.sh  ../bose_2.out'
    start_time = time.time()
    with open('out-{}'.format(magic_str), "w") as outfile:
        result = subprocess.call(cmd, shell=True, stdout=outfile)
    np.savetxt('time-{}'.format(magic_str), [time.time() - start_time])
    benchmark('prop-rho-eq.dat', magic_str, 'bose')
