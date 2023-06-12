import time
import subprocess
import numpy as np
import json
import sys
import itertools
from deom import benchmark, convert, sort_symmetry, complex_2_json, init_qmd, init_qmd_quad

U = 0
Jeff = 0.3

JEXC = [(0.3, 0.3)]
lamd1_l = [2]
hext_l = [0]
lmax_l = [4]
ferr_l = [1e-10]
OMG = 0.5
beta_npfs_l = [(100, 4)]

for (beta, npfs), (Jz, Jp), lmax, lamd1, hext, ferr in itertools.product(
        beta_npfs_l, JEXC, lmax_l, lamd1_l, hext_l, ferr_l):
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
    expn = np.vstack([expn1, expn1, expn2, expn2]).flatten()
    etal = np.vstack([etal1, etal1, etal2, etal2]).flatten()
    etar = np.vstack([etar1, etar1, etar2, etar2]).flatten()
    etaa = np.vstack([etaa1, etaa1, etaa2, etaa2]).flatten()

    nmodmax = np.zeros(len(expn), dtype=complex)
    mode = np.zeros(len(expn), dtype=int)
    for i in range(nmod):
        mode[i * len(etal1):(i + 1) * len(etal1)] = i

    qmdsa = np.zeros((nmod, 4, 4), dtype=complex)
    # start from +/creation, ++--, up-down-up-down
    qmdsa[0, 0, 1] = 1.0
    qmdsa[0, 2, 3] = 1.0
    qmdsa[1, 0, 2] = 1.0
    qmdsa[1, 1, 3] = -1.0
    qmdsa[2, 1, 0] = 1.0
    qmdsa[2, 3, 2] = 1.0
    qmdsa[3, 2, 0] = 1.0
    qmdsa[3, 3, 1] = -1.0

    qmdsc = np.zeros((nmod, 4, 4), dtype=complex)
    # start from -/annihilation
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

    qmds2 = qmds2[:, :, 1:3, 1:3]

    qmds2a = np.zeros(((nmod, nmod, 2, 2)), dtype=complex)
    qmds2a[0, 2] = qmds2[0, 0].copy()
    qmds2a[0, 3] = qmds2[1, 0].copy()
    qmds2a[1, 2] = qmds2[0, 1].copy()
    qmds2a[1, 3] = qmds2[1, 1].copy()
    qmds2a[2, 0] = -qmds2[0, 0].copy()
    qmds2a[2, 1] = -qmds2[0, 1].copy()
    qmds2a[3, 0] = -qmds2[1, 0].copy()
    qmds2a[3, 1] = -qmds2[1, 1].copy()

    qmds2c = np.zeros(((nmod, nmod, 2, 2)), dtype=complex)
    qmds2c[0, 2] = -qmds2[0, 0].copy()
    qmds2c[0, 3] = -qmds2[0, 1].copy()
    qmds2c[1, 2] = -qmds2[1, 0].copy()
    qmds2c[1, 3] = -qmds2[1, 1].copy()
    qmds2c[2, 0] = qmds2[0, 0].copy()
    qmds2c[2, 1] = qmds2[1, 0].copy()
    qmds2c[3, 0] = qmds2[0, 1].copy()
    qmds2c[3, 1] = qmds2[1, 1].copy()

    qmds2b = np.zeros(((nmod, nmod, 2, 2)), dtype=complex)
    qmds2b[0, 0] = -qmds2[0, 0].copy()
    qmds2b[0, 1] = -qmds2[0, 1].copy()
    qmds2b[1, 0] = -qmds2[1, 0].copy()
    qmds2b[1, 1] = -qmds2[1, 1].copy()
    qmds2b[2, 2] = qmds2[0, 0].copy()
    qmds2b[2, 3] = qmds2[1, 0].copy()
    qmds2b[3, 2] = qmds2[0, 1].copy()
    qmds2b[3, 3] = qmds2[1, 1].copy()

    alph1 = 0

    hams = np.zeros((2, 2), dtype=complex)
    for i_s in range(2):
        hams += qmds2[i_s, i_s] * etal1.sum()

    rho0 = np.zeros((2, 2), dtype=complex)
    rho0[0, 0] = 0.5
    rho0[1, 1] = 0.5

    qmdsa = qmdsa[:, 1:3, 1:3]
    qmdsc = qmdsc[:, 1:3, 1:3]

    trun = 0
    lmax = lmax
    nmax = 10000000

    sdip1 = np.zeros((nmod, 2, 2), dtype=complex)
    sdip2 = np.zeros((nmod, 2, 2), dtype=complex)
    sdip3 = np.zeros((nmod, 2, 2), dtype=complex)
    sdip4 = np.zeros((nmod, 2, 2), dtype=complex)

    pdipa = np.zeros((nmod, 2, 2), dtype=complex)
    pdipc = np.zeros((nmod, 2, 2), dtype=complex)
    pdipa[0, :, :] = 0
    pdipa[1, :, :] = 0
    pdipa[2, :, :] = qmds2[0, 0, :, :]
    pdipa[3, :, :] = qmds2[1, 0, :, :]
    pdipc[0, :, :] = qmds2[0, 0, :, :]
    pdipc[1, :, :] = qmds2[1, 0, :, :]
    pdipc[2, :, :] = 0
    pdipc[3, :, :] = 0

    pdipa1 = np.zeros((nmod, 2, 2), dtype=complex)
    pdipc1 = np.zeros((nmod, 2, 2), dtype=complex)
    pdipa1[0, :, :] = qmds2[0, 0, :, :]
    pdipa1[1, :, :] = qmds2[0, 1, :, :]
    pdipa1[2, :, :] = 0
    pdipa1[3, :, :] = 0
    pdipc1[0, :, :] = 0
    pdipc1[1, :, :] = 0
    pdipc1[2, :, :] = 0
    pdipc1[3, :, :] = 0

    bdip = np.ones(np.shape(expn)[0], dtype=complex)

    twsg = np.zeros((npfs, nmod))
    for i in range(npfs):
        for j in range(nmod):
            twsg[i, j] = i + j * npfs

    w_l = []

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
        "alp1": 0,
        "equilibrium": {
            "sc2": False,
            "dt-method": True,
            "OMG": OMG,
            "ti": 0,
            "tf": 10,
            "dt": 0.1,
            "backup": True,
        },
        "read_rho0": True,
        "rho0": complex_2_json(rho0),
        "expn": complex_2_json(expn),
        "ham1": complex_2_json(hams),
        "coef_abs": complex_2_json(etaa),
        "spectrum": True,
        "spectrum-data": {
            "dipole": {
                "sdip_cub": complex_2_json(sdip1),
                "bdip0_cub": complex_2_json(bdip),
                "bdip1_cub": complex_2_json(bdip),
                "bdip2_cub": complex_2_json(bdip),
                "renormalize": complex_2_json(sdip1),
            },
            "dipole1": {
                "sdip_cub": complex_2_json(sdip2),
                "bdip0_cub": complex_2_json(bdip),
                "bdip1_cub": complex_2_json(bdip),
                "bdip2_cub": complex_2_json(bdip),
                "renormalize": complex_2_json(sdip2),
            },
            "if-time": True,
            "time": {
                "Hei": False,
                "ti": 0,
                "tf": 1,
                "dt": 0.01,
                "filter_ferr": ferr,
                "lcr": 'c',
            },
            "file": "prop-pol-1.dat",
        },
    }

    init_qmd(json_init, qmdsa, qmdsc, mode, 2, etaa, etal, etar)
    init_qmd_quad(json_init, qmds2a, qmds2b, qmds2c, mode,
                  2, len(expn), nmod, etaa, etal, etar)
    init_qmd(json_init["spectrum-data"]["dipole"],
             pdipa, pdipc, mode, 2, etaa, etal, etar)
    init_qmd(json_init["spectrum-data"]["dipole1"],
             pdipa1, pdipc1, mode, 2, etaa, etal, etar)

    magic_str = '{}-corr-dt'.format(ferr)
    with open('input.json', 'w') as f:
        json.dump(json_init, f, indent=4, default=convert)
    cmd = r'export OMP_NUM_THREADS={}'.format(
        8 if len(sys.argv) == 1 else sys.argv[1])
    cmd += '&&' + r'/usr/local/bin/jemalloc.sh  ../fermi_quad_2.out'
    start_time = time.time()
    with open('out-{}'.format(magic_str), "w") as outfile:
        result = subprocess.call(cmd, shell=True, stdout=outfile)
    np.savetxt('time-{}'.format(magic_str), [time.time() - start_time])
    benchmark('prop-pol-1.dat', magic_str, 'fermi_quad')
