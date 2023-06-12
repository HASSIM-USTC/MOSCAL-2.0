import time
import numpy as np
import subprocess
import itertools
import json
import sys
from deom import benchmark, convert, complex_2_json, init_qmd, init_qmd_quad


for npsd, lmax, (theta, eta), init, ferr, dt in itertools.product(
    [1],  # npsd
    [100],  # lmax
    [(3 / 4, 0.25)],  # theta, eta
    [1],  # init, 1 for absorption, 0 for emission
    [1e-12],  # ferr
    [0.01],  # dt
):
    nmod = 1
    nsys = 2
    weg = 1
    temp = 1
    OMG = 1
    beta = 1 / temp
    omgb = 1
    ti_gam = 15
    ti_eta = 10
    sdip_val, bdip0_val, bdip1_val, bdip2_val = 1, 1, 1, 1
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
    ],
        dtype=complex)
    etal = np.array([
        0.2703547700830425 - 0.167441733506598j,
        0.8145601320944829 + 0.167441733506598j
    ],
        dtype=complex)
    etar = np.array([
        0.8145601320944829 - 0.167441733506598j,
        0.2703547700830425 + 0.167441733506598j
    ],
        dtype=complex)
    etaa = np.array([0.3180069744932481, 0.831591812680642])

    mode = np.zeros_like(expn, dtype=int)

    hams = np.zeros((nsys, nsys), dtype=complex)
    hams[0, 0] = 1
    hams[0, 1] = 1
    hams[1, 0] = 1

    qmds = np.zeros((nmod, nsys, nsys), dtype=complex)
    if (init == 1):
        qmds[0, 0, 0] = 1
    if (init == 0):
        qmds[0, 1, 1] = 1

    qmd2 = np.zeros((nmod, nmod, 2, 2), dtype=complex)
    for i in range(nmod):
        for j in range(nmod):
            qmd2[i, j, :, :] = qmds[i, :, :]

    nmax = 200000

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
        "filter": True,
        "equilibrium": {
            "sc2": False,
            "dt-method": True,
            "ti": 0,
            "tf": 10,
            "dt": dt,
            "backup": True,
        },
        "expn": complex_2_json(expn),
        "ham1": complex_2_json(hams),
        "coef_abs": complex_2_json(etaa),
        "renormalize": complex_2_json((alp0 + alp2 * np.sum(etal)) * qmds),
    }

    init_qmd(json_init, qmds, qmds, mode, nsys, etaa, etal, etar)
    init_qmd_quad(json_init, qmd2, qmd2, qmd2, mode,
                  nsys, len(expn), nmod, etaa, etal, etar)

    magic_str = '{}-dt'.format(ferr)
    with open('input.json', 'w') as f:
        json.dump(json_init, f, indent=4, default=convert)
    cmd = r'export OMP_NUM_THREADS={}'.format(
        8 if len(sys.argv) == 1 else sys.argv[1])
    cmd += '&&' + r'/usr/local/bin/jemalloc.sh  ../bose_quad_2.out'
    start_time = time.time()
    with open('out-{}'.format(magic_str), "w") as outfile:
        result = subprocess.call(cmd, shell=True, stdout=outfile)
    np.savetxt('time-{}'.format(magic_str), [time.time() - start_time])
    benchmark('prop-rho-eq.dat', magic_str, 'bose_quad')
