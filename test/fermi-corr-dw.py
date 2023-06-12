import time
import subprocess
import numpy as np
import json
import itertools
import sys
from deom import benchmark, convert, sort_symmetry, complex_2_json, init_qmd


beta_l = [100]
lamd1_u_l = [(0.1, 1.2)]
lmax_l = [4]
ferr_l = [1e-8]
OMG = 0.5
npfs = 5
iterlist = itertools.product(beta_l, lmax_l, lamd1_u_l, ferr_l)

for beta, lmax, (lamd1, U), ferr in iterlist:
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

    nmodmax = np.zeros(len(expn), dtype=float)
    mode = np.zeros(len(expn), dtype=int)
    for i in range(nmod):
        mode[i * len(etal1):(i + 1) * len(etal1)] = i

    hams = np.zeros((4, 4), dtype=complex)
    hams[0, 0] = 0
    hams[1, 1] = -U / 2
    hams[2, 2] = -U / 2
    hams[3, 3] = 0

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

    lmax = lmax
    nmax = 10000000

    sdip1 = np.zeros((nmod, 4, 4), dtype=complex)
    sdip2 = np.zeros((nmod, 4, 4), dtype=complex)
    sdip3 = np.zeros((nmod, 4, 4), dtype=complex)
    sdip4 = np.zeros((nmod, 4, 4), dtype=complex)

    sdip1[0, :, :] = qmdsc[0, :, :]
    sdip2[0, :, :] = qmdsa[0, :, :]
    sdip3[0, :, :] = qmdsc[0, :, :]
    sdip4[0, :, :] = qmdsa[0, :, :]

    pdipa = np.zeros((nmod, 4, 4), dtype=complex)
    pdipc = np.zeros((nmod, 4, 4), dtype=complex)
    bdip = np.ones(np.shape(expn)[0], dtype=complex)

    twsg = np.zeros((npfs, nmod))
    for i in range(npfs):
        for j in range(nmod):
            twsg[i, j] = i + j * npfs

    w_l = np.linspace(-0.1, 0, 3)

    json_init = {
        "syl": {
            "OMG": OMG,
            "nind": npfs,
            "lwsg": nmod,
            "twsg": list(twsg.flatten())
        },
        "nmax": nmax,
        "lmax": lmax,
        "ferr": ferr,
        "filter": False,
        "nind": len(expn),
        "nmod": nmod,
        "equilibrium": {
            "sc2": True,
            "dt-method": False,
            "OMG": OMG,
            "ti": 0,
            "tf": 50,
            "dt": 0.1,
            "backup": True,
        },
        "expn": complex_2_json(expn),
        "ham1": complex_2_json(hams),
        "coef_abs": complex_2_json(etaa),
        "spectrum": True,
        "spectrum-data": {
            "dipole": {
                "sdip_cub": complex_2_json(sdip1),
                "bdip1_cub": complex_2_json(bdip),
            },
            "dipole1": {
                "sdip_cub": complex_2_json(sdip2),
                "bdip1_cub": complex_2_json(bdip),
            },
            "if-time": False,
            "frequency": {
                "lcr": 'c',
                "Hei": False,
                "noise": False,
                "OMG": OMG,
                "w_len": len(w_l),
                "w_l": list(w_l),
                "ferr": 1e-8,
                "filter_ferr": ferr,
                "step": 10,
            },
            "file": "prop-pol-1.dat",
        },
    }

    init_qmd(json_init, qmdsa, qmdsc, mode, 4, etaa, etal, etar)
    init_qmd(json_init["spectrum-data"]["dipole"],
             pdipa, pdipc, mode, 4, etaa, etal, etar)
    init_qmd(json_init["spectrum-data"]["dipole1"],
             pdipa, pdipc, mode, 4, etaa, etal, etar)

    magic_str = '{}-corr-dw'.format(ferr)
    with open('input.json', 'w') as f:
        json.dump(json_init, f, indent=4, default=convert)
    cmd = r'export OMP_NUM_THREADS={}'.format(
        8 if len(sys.argv) == 1 else sys.argv[1])
    cmd += '&&' + r'/usr/local/bin/jemalloc.sh  ../fermi_4.out'
    start_time = time.time()
    with open('out-{}'.format(magic_str), "w") as outfile:
        result = subprocess.call(cmd, shell=True, stdout=outfile)
    np.savetxt('time-{}'.format(magic_str), [time.time() - start_time])
    benchmark('prop-pol-1.dat', magic_str, 'fermi', if_np=False)
