import time
import subprocess
import numpy as np
import json
import itertools
import sys
from deom import benchmark, convert, sort_symmetry, complex_2_json, init_qmd, init_qmd_quad


beta_l = [100]
JEXC = [(2, 2)]
lamd1_u_l = [(0.5, 0.8)]
lmax_l = [4]
ferr_l = [1e-10]
OMG = 0.5
npfs_l = [5]

for beta, (Jz, Jp), lmax, (lamd1, U), ferr, npfs in itertools.product(
        beta_l, JEXC, lmax_l, lamd1_u_l, ferr_l, npfs_l):
    file_path = "PFS/fermi-1-1/f-{}-comp/".format(beta)
    file_path += "{}" + "_{}".format(npfs - 1)
    etal1 = lamd1 * np.loadtxt(file_path.format("etal1"), dtype=complex)
    expn1 = np.loadtxt(file_path.format("expn1"), dtype=complex)
    etal2 = lamd1 * np.loadtxt(file_path.format("etal2"), dtype=complex)
    expn2 = np.loadtxt(file_path.format("expn2"), dtype=complex)

    expn1, etal1, etar1, etaa1 = sort_symmetry(etal1, expn1, False)
    expn2, etal2, etar2, etaa2 = sort_symmetry(etal2, expn2, False)

    dt = 0.01
    nmod = 4
    alp1 = 0
    lmax = lmax
    nmax = 10000000
    expn = np.vstack([expn1, expn1, expn2, expn2]).flatten()
    etal = np.vstack([etal1, etal1, etal2, etal2]).flatten()
    etar = np.vstack([etar1, etar1, etar2, etar2]).flatten()
    etaa = np.vstack([etaa1, etaa1, etaa2, etaa2]).flatten()

    nmodmax = np.zeros(len(expn), dtype=float)
    mode = np.zeros(len(expn), dtype=int)
    for i in range(nmod):
        mode[i * len(etal1):(i + 1) * len(etal1)] = i

    hams0 = np.zeros((4, 4), dtype=complex)

    qmdsa = np.zeros((nmod, 4, 4), dtype=complex)
    qmdsa[0, 0, 1] = 1.0
    qmdsa[0, 2, 3] = 1.0
    qmdsa[1, 0, 2] = 1.0
    qmdsa[1, 1, 3] = -1.0
    qmdsa[2, 1, 0] = 1.0
    qmdsa[2, 3, 2] = 1.0
    qmdsa[3, 2, 0] = 1.0
    qmdsa[3, 3, 1] = -1.0

    qmdsc = np.zeros((nmod, 4, 4), dtype=complex)
    qmdsc[0, 1, 0] = 1.0
    qmdsc[0, 3, 2] = 1.0
    qmdsc[1, 2, 0] = 1.0
    qmdsc[1, 3, 1] = -1.0
    qmdsc[2, 0, 1] = 1.0
    qmdsc[2, 2, 3] = 1.0
    qmdsc[3, 0, 2] = 1.0
    qmdsc[3, 1, 3] = -1.0

    qmds2 = np.zeros((2, 2, 4, 4), dtype=complex)

    qmds2[0, 0] = Jz * (qmdsa[0].dot(qmdsa[2]) - qmdsa[1].dot(qmdsa[3]))
    qmds2[1, 1] = -qmds2[0, 0].copy()
    qmds2[0, 1] = 2 * qmdsa[1].dot(qmdsa[2]) * Jp
    qmds2[1, 0] = 2 * qmdsa[0].dot(qmdsa[3]) * Jp

    qmds2 = qmds2 / 4

    qmds2a = np.zeros(((nmod, nmod, 4, 4)), dtype=complex)
    qmds2a[0, 2] = qmds2[0, 0].copy()
    qmds2a[0, 3] = qmds2[1, 0].copy()
    qmds2a[1, 2] = qmds2[0, 1].copy()
    qmds2a[1, 3] = qmds2[1, 1].copy()
    qmds2a[2, 0] = -qmds2[0, 0].copy()
    qmds2a[2, 1] = -qmds2[0, 1].copy()
    qmds2a[3, 0] = -qmds2[1, 0].copy()
    qmds2a[3, 1] = -qmds2[1, 1].copy()

    qmds2c = np.zeros(((nmod, nmod, 4, 4)), dtype=complex)
    qmds2c[0, 2] = -qmds2[0, 0].copy()
    qmds2c[0, 3] = -qmds2[0, 1].copy()
    qmds2c[1, 2] = -qmds2[1, 0].copy()
    qmds2c[1, 3] = -qmds2[1, 1].copy()
    qmds2c[2, 0] = qmds2[0, 0].copy()
    qmds2c[2, 1] = qmds2[1, 0].copy()
    qmds2c[3, 0] = qmds2[0, 1].copy()
    qmds2c[3, 1] = qmds2[1, 1].copy()

    qmds2b = np.zeros(((nmod, nmod, 4, 4)), dtype=complex)
    qmds2b[0, 0] = -qmds2[0, 0].copy()
    qmds2b[0, 1] = -qmds2[0, 1].copy()
    qmds2b[1, 0] = -qmds2[1, 0].copy()
    qmds2b[1, 1] = -qmds2[1, 1].copy()
    qmds2b[2, 2] = qmds2[0, 0].copy()
    qmds2b[2, 3] = qmds2[1, 0].copy()
    qmds2b[3, 2] = qmds2[0, 1].copy()
    qmds2b[3, 3] = qmds2[1, 1].copy()

    hams = np.zeros((4, 4), dtype=complex)
    hams = hams0.copy()
    for i_s in range(2):
        hams += qmds2[i_s, i_s] * etal1.sum()

    twsg = np.zeros((npfs, nmod))
    for i in range(npfs):
        for j in range(nmod):
            twsg[i, j] = i + j * npfs

    json_init = {
        "syl": {
            "OMG": OMG,
            "nind": npfs,
            "lwsg": nmod,
            "twsg": list(twsg.flatten())
        },
        "nmax": nmax,
        "lmax": lmax,
        "filter": False,
        "nind": len(expn),
        "nmod": nmod,
        "inistate": 2,
        "alp1": alp1,
        "equilibrium": {
            "sc2": False,
            "dt-method": True,
            "OMG": OMG,
            "ti": 0,
            "tf": 10,
            "dt": 0.1,
            "backup": True,
        },
        "expn": complex_2_json(expn),
        "ham1": complex_2_json(hams),
        "coef_abs": complex_2_json(etaa),
    }

    init_qmd(json_init, qmdsa, qmdsc, mode, 4, etaa, etal, etar)
    init_qmd_quad(json_init, qmds2a, qmds2b, qmds2c, mode,
                  4, len(expn), nmod, etaa, etal, etar)

    magic_str = '{}-fermi-quad-dt'.format(ferr)
    with open('input.json', 'w') as f:
        json.dump(json_init, f, indent=4, default=convert)
    cmd = r'export OMP_NUM_THREADS={}'.format(
        8 if len(sys.argv) == 1 else sys.argv[1])
    cmd += '&&' + r'/usr/local/bin/jemalloc.sh  ../fermi_quad_4.out'
    start_time = time.time()
    with open('out-{}'.format(magic_str), "w") as outfile:
        result = subprocess.call(cmd, shell=True, stdout=outfile)
    np.savetxt('time-{}'.format(magic_str), [time.time() - start_time])
    benchmark('prop-rho-eq.dat', magic_str, 'fermi_quad')
