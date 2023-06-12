import time
import subprocess
import numpy as np
import json
import itertools
import pandas as pd
from termcolor import colored
from aux_function import fastread_np, sum_exptontial_spectrum, convert, fermi_ado_number

beta_l = [100]
JEXC = [(2, 2)]
lamd1_u_l = [(0.5, 0.8)]
lmax_l = [4]
ferr_l = [1e-10]
OMG = 0.25
npfs_l = [7]

for beta, (Jz, Jp), lmax, (lamd1, U), ferr, npfs in itertools.product(
        beta_l, JEXC, lmax_l, lamd1_u_l, ferr_l, npfs_l):
    etal1 = lamd1 * np.loadtxt("PFS/fermi-1-1/f-{}-comp/etal1_{}".format(
        beta, npfs - 1),
                               dtype=complex)
    expn1 = np.loadtxt("PFS/fermi-1-1/f-{}-comp/expn1_{}".format(
        beta, npfs - 1),
                       dtype=complex)
    etal2 = lamd1 * np.loadtxt("PFS/fermi-1-1/f-{}-comp/etal2_{}".format(
        beta, npfs - 1),
                               dtype=complex)
    expn2 = np.loadtxt("PFS/fermi-1-1/f-{}-comp/expn2_{}".format(
        beta, npfs - 1),
                       dtype=complex)

    expn_imag_sort = np.argsort(np.abs(np.imag(expn1)))[::-1]
    expn_imag = np.sort(np.abs(np.imag(expn1)))[::-1]
    expn1 = expn1[expn_imag_sort]
    etal1 = etal1[expn_imag_sort]
    etar1 = etal1[expn_imag_sort]
    expn_val_cc = np.where(expn1[expn_imag > 1e-14])[0]
    etaa1 = np.zeros(len(etal1), dtype=float)
    for ii in range(0, len(expn_val_cc), 2):
        even_i = ii
        odd_i = ii + 1
        etar1[even_i] = np.conj(etal1[odd_i])
        etar1[odd_i] = np.conj(etal1[even_i])
        etaa1[even_i] = np.abs(etal1[even_i])
        etaa1[odd_i] = np.abs(etal1[odd_i])
    for ii in range(len(expn_val_cc), npfs):
        even_i = ii
        etar1[even_i] = np.conj(etal1[even_i])
        etaa1[even_i] = np.abs(etal1[even_i])

    expn_imag_sort = np.argsort(np.abs(np.imag(expn2)))[::-1]
    expn_imag = np.sort(np.abs(np.imag(expn1)))[::-1]
    expn2 = expn2[expn_imag_sort]
    etal2 = etal2[expn_imag_sort]
    etar2 = etal2[expn_imag_sort]
    expn_val_cc = np.where(expn2[expn_imag > 1e-14])[0]
    etaa2 = np.zeros(len(etal2), dtype=float)
    for ii in range(0, len(expn_val_cc), 2):
        even_i = ii
        odd_i = ii + 1
        etar2[even_i] = np.conj(etal2[odd_i])
        etar2[odd_i] = np.conj(etal2[even_i])
        etaa2[even_i] = np.abs(etal2[even_i])
        etaa2[odd_i] = np.abs(etal2[odd_i])
    for ii in range(len(expn_val_cc), npfs):
        even_i = ii
        etar2[even_i] = np.conj(etal2[even_i])
        etaa2[even_i] = np.abs(etal2[even_i])

    # len_ = 1000000
    # spe_wid = 20
    # w = np.linspace(-spe_wid, spe_wid, len_)
    # res_J1 = np.zeros(len(w), dtype=complex)
    # fit_J(w, res_J1, expn1, etal1, 1)
    # plt.plot(w, (res_J1.real), "g", label="PFS")
    # res_J1 = np.zeros(len(w), dtype=complex)
    # fit_J(w, res_J1, expn2, etal2, -1)
    # plt.plot(w, (res_J1.real), "r", label="PFS")
    # plt.savefig("spe.pdf")

    dt = 0.01
    nmod = 4
    expn = np.append(np.append(np.append(expn1, expn1), expn2), expn2)
    etal = np.append(np.append(np.append(etal1, etal1), etal2), etal2)
    etar = np.append(np.append(np.append(etar1, etar1), etar2), etar2)
    etaa = np.append(np.append(np.append(etaa1, etaa1), etaa2), etaa2)

    nmodmax = np.zeros(len(expn), dtype=float)
    mode = np.zeros(len(expn), dtype=float)
    mode[:len(etal1)] = 0
    mode[len(etal1):2 * len(etal1)] = 1
    mode[2 * len(etal1):3 * len(etal1)] = 2
    mode[3 * len(etal1):] = 3

    hams0 = np.zeros((4, 4), dtype=complex)
    # hams0[0, 0] = 0
    # hams0[1, 1] = -U / 2
    # hams0[2, 2] = -U / 2
    # hams0[3, 3] = 0

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

    alp1 = 0

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

    w_l = np.linspace(-0.3, 0, 31)
    # w_l = []
    # w_l = np.linspace(-0.2, 0, 81)[1:]
    # w_l = np.append(np.linspace(-1, -0.1, 41), np.linspace(-0.1, 0, 21)[1:])

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
        "ferr": ferr,
        # "filter": True,
        "filter": False,
        "nind": len(expn),
        "nmod": nmod,
        "nmodmax": 1,
        "inistate": 2,
        "alp1": alp1,
        "equilibrium": {
            "sc2": True,
            "cgsb": False,
            "dt-method": False,
            "re_tree": 10000,
            "filter": 10000,
            "OMG": OMG,
            "ti": 0,
            "tf": 100,
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
            "lcr": 'c',
            "Hei": False,
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
    cmd = r'mv prop-pol1-1.dat prop-pol1-1.dat-{}-fermi-quad-sc2-dt'.format(
        ferr)
    print(time.time() - start_time)
    result = subprocess.call(cmd, shell=True)

    data1 = fastread_np(
        './result/prop-pol1-1.dat-{}-fermi-quad-sc2-dt'.format(ferr))
    data2 = fastread_np('./prop-pol1-1.dat-{}-fermi-quad-sc2-dt'.format(ferr))
    result = (np.sum(np.abs(data1 - data2)))
    print(
        result,
        colored('PASSED', 'green') if float(result) < 1e-4 else colored(
            'FAILED', 'red'))
