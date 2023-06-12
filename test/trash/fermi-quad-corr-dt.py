import time
import subprocess
import numpy as np
import json
import itertools
from termcolor import colored
from aux_function import fastread_np, convert, sort_symmetry

beta_l = [4]
JEXC = [(2, 2)]
lamd1_u_l = [(0.5, 0.8)]
lmax_l = [4]
ferr_l = [1e-20]
OMG = 0.5
npfs_l = [2]

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

    print(qmds2a)
    print(qmds2b)
    print(qmds2c)

    alph1 = 0

    hams = np.zeros((4, 4), dtype=complex)
    hams = hams0.copy()
    for i_s in range(2):
        hams += qmds2[i_s, i_s] * etal1.sum()

    # hidx
    trun = 0
    lmax = lmax
    nmax = 10000000

    sdip1 = np.zeros((nmod, 4, 4), dtype=float)
    sdip2 = np.zeros((nmod, 4, 4), dtype=float)
    sdip3 = np.zeros((nmod, 4, 4), dtype=float)
    sdip4 = np.zeros((nmod, 4, 4), dtype=float)

    sdip1[0, :, :] = qmdsc[0, :, :]
    sdip2[0, :, :] = qmdsa[0, :, :]
    sdip3[0, :, :] = qmdsc[0, :, :]
    sdip4[0, :, :] = qmdsa[0, :, :]

    pdipa = np.zeros((nmod, 4, 4), dtype=float)
    pdipc = np.zeros((nmod, 4, 4), dtype=float)
    bdip = np.ones(np.shape(expn)[0], dtype=float)

    twsg = np.zeros((npfs, nmod))
    for i in range(npfs):
        for j in range(nmod):
            twsg[i, j] = i + j * npfs

    twsa = np.zeros((nmod // 2, npfs))
    tsad = np.zeros((nmod // 2, npfs))
    for i in range(nmod // 2):
        for j in range(npfs):
            twsa[i, j] = i * npfs + j
            tsad[i, j] = i * npfs + j + npfs * (nmod // 2)

    w_l = []

    amax = 2 * np.ones(nmod)
    dmax = 2 * np.ones(nmod)

    jsonInit = {
        "syl": {
            "OMG": OMG,
            "nind": npfs,
            "lwsg": nmod,
            "twsg": list(twsg.flatten())
        },
        "sym": {
            "nind": nmod // 2,
            "lwsa": npfs,
            "lsad": npfs,
            "twsa": list(twsa.flatten()),
            "tsad": list(tsad.flatten()),
            "amax": list(amax.flatten()),
            "dmax": list(dmax.flatten()),
        },
        "nmax": nmax,
        "lmax": lmax,
        "nsys": 4,
        "nham": 4,
        "filter": False,
        "nind": len(expn),
        "nmod": nmod,
        "nmodmax": 1,
        "inistate": 2,
        "equilibrium": {
            "sc2": False,
            "dt-method": True,
            "re_tree": 10000,
            "filter": 10000,
            "OMG": OMG,
            "ti": 0,
            "tf": 1000,
            "dt": 0.1,
            "backup": True,
        },
        "cal": {
            "corr": True,
            "11": True,
            "22": False,
            "21": False,
            "12": False,
        },
        "corr": {
            "der": False,
            "noise": False,
            "ti": 0,
            "tf": 1000,
            "dt": 0.1,
            "filter_ferr": ferr,
            "lcr": 'c',
        },
        "corr-w": {
            "Hei": True,
            "noise": False,
            "OMG": OMG,
            "w_len": len(w_l),
            "w_l": list(w_l),
            "ferr": 1e-3,
            "filter_ferr": ferr,
            "step": 10,
            "re_tree": 10000,
            "filter": 10000,
            "ti": 0,
            "tf": 10,
            "dt": 0.1,
        },
        "expn": {
            "real": list(np.real(expn.flatten())),
            "imag": list(np.imag(expn.flatten()))
        },
        "ham1": {
            "real": list(np.real(hams.flatten())),
            "imag": list(np.imag(hams.flatten()))
        },
        "alph1": alph1,
        "qmd1a": {
            "real": list(np.real(qmdsa.flatten())),
            "imag": list(np.imag(qmdsa.flatten()))
        },
        "qmd1c": {
            "real": list(np.real(qmdsc.flatten())),
            "imag": list(np.imag(qmdsc.flatten()))
        },
        "qmd2a": {
            "real": list(np.real(qmds2a.flatten())),
            "imag": list(np.imag(qmds2a.flatten()))
        },
        "qmd2b": {
            "real": list(np.real(qmds2b.flatten())),
            "imag": list(np.imag(qmds2b.flatten()))
        },
        "qmd2c": {
            "real": list(np.real(qmds2c.flatten())),
            "imag": list(np.imag(qmds2c.flatten()))
        },
        "modLabel": {
            "real": list(np.real(mode.flatten())),
            "imag": list(np.imag(mode.flatten()))
        },
        "nmodmaxLabel": {
            "real": list(np.real(nmodmax.flatten())),
            "imag": list(np.imag(nmodmax.flatten()))
        },
        "delt_res": {
            "real": list(np.zeros(nmod)),
            "imag": list(np.zeros(nmod))
        },
        "coef_abs": {
            "real": list(np.real(etaa.flatten())),
            "imag": list(np.imag(etaa.flatten()))
        },
        "coef_lft": {
            "real": list(np.real(etal.flatten())),
            "imag": list(np.imag(etal.flatten()))
        },
        "coef_rht": {
            "real": list(np.real(etar.flatten())),
            "imag": list(np.imag(etar.flatten()))
        },
        "dipole": {
            "sdip_cub": {
                "real": list(np.real(sdip1.flatten())),
                "imag": list(np.imag(sdip1.flatten()))
            },
            "pdipa_cub": {
                "real": list(np.real(pdipa.flatten())),
                "imag": list(np.imag(pdipa.flatten()))
            },
            "pdipc_cub": {
                "real": list(np.real(pdipc.flatten())),
                "imag": list(np.imag(pdipc.flatten()))
            },
            "bdip_cub": {
                "real": list(np.real(bdip.flatten())),
                "imag": list(np.imag(bdip.flatten()))
            }
        },
        "dipole1": {
            "sdip_cub": {
                "real": list(np.real(sdip2.flatten())),
                "imag": list(np.imag(sdip2.flatten()))
            },
            "pdipa_cub": {
                "real": list(np.real(pdipa.flatten())),
                "imag": list(np.imag(pdipa.flatten()))
            },
            "pdipc_cub": {
                "real": list(np.real(pdipc.flatten())),
                "imag": list(np.imag(pdipc.flatten()))
            },
            "bdip_cub": {
                "real": list(np.real(bdip.flatten())),
                "imag": list(np.imag(bdip.flatten()))
            }
        },
        "dipole2": {
            "sdip_cub": {
                "real": list(np.real(sdip3.flatten())),
                "imag": list(np.imag(sdip3.flatten()))
            },
            "pdipa_cub": {
                "real": list(np.real(pdipa.flatten())),
                "imag": list(np.imag(pdipa.flatten()))
            },
            "pdipc_cub": {
                "real": list(np.real(pdipc.flatten())),
                "imag": list(np.imag(pdipc.flatten()))
            },
            "bdip_cub": {
                "real": list(np.real(bdip.flatten())),
                "imag": list(np.imag(bdip.flatten()))
            }
        },
        "dipole3": {
            "sdip_cub": {
                "real": list(np.real(sdip4.flatten())),
                "imag": list(np.imag(sdip4.flatten()))
            },
            "pdipa_cub": {
                "real": list(np.real(pdipa.flatten())),
                "imag": list(np.imag(pdipa.flatten()))
            },
            "pdipc_cub": {
                "real": list(np.real(pdipc.flatten())),
                "imag": list(np.imag(pdipc.flatten()))
            },
            "bdip_cub": {
                "real": list(np.real(bdip.flatten())),
                "imag": list(np.imag(bdip.flatten()))
            }
        }
    }

    with open('input.json', 'w') as f:
        json.dump(jsonInit, f, indent=4, default=convert)
    cmd = r'export OMP_NUM_THREADS=8' + '&&'
    cmd += r'/usr/local/bin/jemalloc.sh ../fermi_quad_4.out'
    start_time = time.time()
    with open('out-{}-{}-{}-dt'.format(Jz, Jp, npfs), "w") as outfile:
        result = subprocess.call(cmd, shell=True, stdout=outfile)
    np.savetxt('time-{}-{}-{}-dt'.format(Jz, Jp, npfs),
               [time.time() - start_time])
    print(time.time() - start_time)
    cmd = r'mv prop-pol1-1.dat prop-pol1-1.dat-{}-fermi-quad-dt'.format(ferr)
    result = subprocess.call(cmd, shell=True)

    data1 = fastread_np(
        './result/prop-pol1-1.dat-{}-fermi-quad-dt'.format(ferr))
    data2 = fastread_np('./prop-pol1-1.dat-{}-fermi-quad-dt'.format(ferr))
    result = (np.sum(np.abs(data1 - data2)))
    print(
        result,
        colored('PASSED', 'green') if float(result) < 1e-4 else colored(
            'FAILED', 'red'))
