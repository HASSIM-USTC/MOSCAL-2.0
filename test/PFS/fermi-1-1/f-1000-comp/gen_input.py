import math
from numpy import linalg as LA
import re
import os
import time
from io import StringIO
import pandas as pd
import numpy as np
import sympy as sp
import json
import gc
from functools import partial
import numpy as np
from matplotlib import rc
import matplotlib.pyplot as plt
import itertools


def convert(o):
    if isinstance(o, np.int64): return int(o)
    raise TypeError


def fFeim(x, pole, resi, sigma=1):
    return 1 / 2 - sigma * sum(2.0 * resi[i] * x / (x**2 + pole[i]**2)
                               for i in range(len(pole)))


def tseig(D, E):
    mat = np.diag(E, -1) + np.diag(D, 0) + np.diag(E, 1)
    return -np.sort(-np.linalg.eigvalsh(mat))


def MSD(N, BoseFermi=1):
    if BoseFermi == 1:
        pole = np.array([2 * (i + 1) * np.pi for i in range(N)])
        resi = np.ones(N, dtype=float)
        return pole, resi
    elif BoseFermi == 2:
        pole = np.array([(2 * i + 1) * np.pi for i in range(N)])
        resi = np.ones(N, dtype=float)
        return pole, resi


def PSD(N, BoseFermi=1, pade=0):
    if N < 0 or BoseFermi < 1 or BoseFermi > 2 or pade < 0 or pade > 3:
        raise ValueError("N or BoseFermi or pade has wrong value!")

    if pade == 0:
        return MSD(N, BoseFermi)
    elif pade == 1 or pade == 2:
        pole, resi = [], []
        if N > 0:
            M = 2 * N + pade // 2
            temp = 3.0 if BoseFermi == 1 else 1.0
            diag = np.zeros(M, dtype=float)
            doff = np.array([
                1.0 / math.sqrt((temp + 2.0 * i) * (temp + 2.0 * (i + 1)))
                for i in range(M - 1)
            ])
            pole = 2.0 / tseig(diag, doff)[:N]
            pol2 = np.array([x * x for x in pole])
            M -= 1
            temp = 5.0 if BoseFermi == 1 else 3.0
            diag = np.zeros(M, dtype=float)
            doff = np.array([
                1.0 / math.sqrt((temp + 2.0 * i) * (temp + 2.0 * (i + 1)))
                for i in range(M - 1)
            ])
            M //= 2
            eig2 = np.power(2.0 / tseig(diag, doff)[:M], 2)
            scaling = 0.0
            if BoseFermi == 1:
                scaling = N*(2.0*N+3.0) if pade == 1 else 1.0 / \
                    (4.0*(N+1.0)*(2.0*N+3.0))
            elif BoseFermi == 2:
                scaling = N*(2.0*N+1.0) if pade == 1 else 1.0 / \
                    (4.0*(N+1.0)*(2.0*N+1.0))
            resi = np.zeros(N, dtype=float)
            for j in range(N):
                if pade == 2:
                    temp = 0.5 * scaling * (eig2[j] - pol2[j])
                elif pade == 1:
                    if j == N - 1:
                        temp = 0.5 * scaling
                    else:
                        temp = 0.5*scaling * \
                            (eig2[j]-pol2[j])/(pol2[N-1]-pol2[j])
                for k in range(M):
                    temp *= (eig2[k]-pol2[j]) / \
                        (pol2[k]-pol2[j]) if k != j else 1.0
                resi[j] = temp
        rn, tn = 0.0, 0.0
        if BoseFermi == 1 and pade == 2:
            rn = 1.0 / (4.0 * (N + 1.0) * (2.0 * N + 3.0))
        return pole, resi
    elif pade == 3:
        Np1 = N + 1
        temp = 3.0 if BoseFermi == 1 else 1.0
        d = np.empty(2 * Np1, dtype=float)
        d[0] = 0.25 / temp
        d[-1] = -4.0 * (N + 1.0) * (N + 1.0) * (temp + 2 * N) * (
            temp + 2 * N) * (temp + 4 * N + 2.0)
        for i in range(1, Np1):
            d[2*i-1] = -4.0*i*i*(temp+2.0*i-2.0) * \
                (temp+2.0*i-2.0)*(temp+4.0*i-2.0)
            d[2 * i] = -0.25 * (temp + 4.0 * i) / i / (i + 1) / (
                temp + 2.0 * i - 2.0) / (temp + 2.0 * i)
        sumd2 = np.empty(Np1, dtype=float)
        sumd2[0] = d[1]
        for i in range(1, Np1):
            sumd2[i] = sumd2[i - 1] + d[2 * i + 1]
        tn = 0.25 / sumd2[-1]
        rn = sum(d[2 * i] * (4.0 * tn *
                             (sumd2[-1] - sumd2[i - 1]))**2 if i > 0 else d[2 *
                                                                            i]
                 for i in range(Np1))
        M = 2 * N + 1
        diag = np.zeros(M, dtype=float)
        doff = np.array(
            [1.0 / math.sqrt(d[i + 1] * d[i + 2]) for i in range(M - 1)])
        pole = 2.0 / tseig(diag, doff)[:N]
        resi = np.zeros(N, dtype=float)
        for j in range(N):
            scaling = pole[j] * pole[j]
            r0, t1 = 0.0, 0.25 / d[1]
            eta0, eta1, eta2 = 0.0, 0.5, 0.0
            for i in range(Np1):
                r1 = t1 if (i == j
                            or i == N) else t1 / (pole[i] * pole[i] - scaling)
                r2 = 2.0*math.sqrt(abs(r1)) if r1 > 0 else - \
                    2.0*math.sqrt(abs(r1))
                r1 = 2.0 * math.sqrt(abs(r1))
                eta2 = d[2 * i] * r1 * eta1 - 0.25 * r1 * r0 * scaling * eta0
                eta0 = eta1
                eta1 = eta2
                eta2 = d[2 * i +
                         1] * r2 * eta1 - 0.25 * r2 * r1 * scaling * eta0
                eta0 = eta1
                eta1 = eta2
                r0 = r2
                if i != N:
                    t1 = sumd2[i] / sumd2[i + 1]
            resi[j] = eta2
        return pole, resi


# in this script, we can decompose any given spectrum, but the sympy format is must been given
# do u like haskell?
# sympy[spe(def by sympy)], dict[sp_para_dict], dict[para_dict], dict[npsd],
# dict[pade] >> np.array[etal], np.array[etar],np.array[etaa], np.array[expn]
def decompose_spe(spe,
                  sp_para_dict,
                  para_dict,
                  condition_dict,
                  npsd,
                  sigma,
                  pade=1):
    numer, denom = sp.cancel(sp.factor(
        spe.subs(condition_dict))).as_numer_denom()
    numer_get_para = (sp.factor(numer)).subs(sp_para_dict)
    denom_get_para = (sp.factor(denom)).subs(sp_para_dict)
    poles = sp.nroots(denom_get_para)
    float(sp.re(poles[0]))

    expn = []
    poles_allplane = np.array([])
    for i in poles:
        i = complex(i)
        if i.imag * sigma > 0:
            expn.append(-sigma * i * 1.J)
        poles_allplane = np.append(poles_allplane, i)

    etal = []
    etar = []
    etaa = []

    expn = np.array(expn)
    expn_imag_sort = np.argsort(np.abs(np.imag(expn)))[::-1]
    expn_val_n_cc = expn

    expn = list(expn[expn_imag_sort])
    pole, resi = PSD(npsd, 2, 1)
    beta = para_dict['beta']
    temp = 1 / beta

    for ii in range(len(expn_val_n_cc)):
        etal.append(
            complex(
                sp.N((sigma * 2.j * numer_get_para / np.multiply.reduce(
                    w_sp - poles_allplane[np.abs(poles_allplane - sigma * 1.J *
                                                 expn_val_n_cc[ii]) > 1e-14])
                      ).subs({w_sp: sigma * 1.j * expn_val_n_cc[ii]}) *
                     fFeim(sigma * 1.J * expn_val_n_cc[ii] / temp, pole, resi,
                           sigma))))
        etar.append(np.conj(etal[-1]))
        etaa.append(np.sqrt(np.abs(etal[-1]) * np.abs(etar[-1])))

    f = numer_get_para / np.multiply.reduce(w_sp - poles_allplane)
    f = sp.lambdify(w_sp, f)

    for inma in range(len(pole)):
        print(pole[inma])
        zomg = sigma * 1.J * pole[inma] * temp
        jsum = np.sum(f(zomg))
        expn.append(pole[inma] * temp)
        etal.append(-2.J * resi[inma] * temp * jsum)
        etar.append(np.conj(etal[-1]))
        etaa.append(np.abs(etal[-1]))

    etal = np.array(etal)
    etar = np.array(etar)
    etaa = np.array(etaa)
    expn = np.array(expn)
    return etal, etar, etaa, expn


def decompose_spe_real(spe,
                       sp_para_dict,
                       para_dict,
                       condition_dict,
                       sigma,
                       pade=1):
    numer, denom = sp.cancel(sp.factor(
        spe.subs(condition_dict))).as_numer_denom()
    numer_get_para = (sp.factor(numer)).subs(sp_para_dict)
    denom_get_para = (sp.factor(denom)).subs(sp_para_dict)
    poles = sp.nroots(denom_get_para)
    float(sp.re(poles[0]))
    print(poles)

    expn = []
    poles_allplane = np.array([])
    for i in poles:
        i = complex(i)
        if i.imag * sigma > 0:
            expn.append(-sigma * i * 1.J)
        poles_allplane = np.append(poles_allplane, i)

    etal = []
    etar = []
    etaa = []
    expn = np.array(expn)
    expn_imag_sort = np.argsort(np.abs(np.imag(expn)))[::-1]
    expn_val_n_cc = expn
    expn = list(expn[expn_imag_sort])

    for ii in range(len(expn_val_n_cc)):
        etal.append(
            complex(
                sp.N((sigma * 1.j * numer_get_para /
                      np.multiply.reduce(w_sp - poles_allplane[np.abs(
                          poles_allplane -
                          sigma * 1.J * expn_val_n_cc[ii]) > 1e-14])).subs(
                              {w_sp: sigma * 1.j * expn_val_n_cc[ii]}))))
        etar.append(np.conj(etal[-1]))
        etaa.append(np.sqrt(np.abs(etal[-1]) * np.abs(etar[-1])))

    etal = np.array(etal)
    etar = np.array(etar)
    etaa = np.array(etaa)
    expn = np.array(expn)
    return etal, etar, etaa, expn


def fit_J(w, res, expn, etal, sigma):
    for i in range(len(etal)):
        res += etal[i] / (expn[i] + sigma * 1.j * w)


def fit_t(t, res, expn, etal):
    for i in range(len(etal)):
        res += etal[i] * np.exp(-expn[i] * t)
    return res


def fit_t_d1(t, res, expn, etal):
    for i in range(len(etal)):
        res += -expn[i] * etal[i] * np.exp(-expn[i] * t)
    return res


def get_Q_h(n, expn_pade1, etal_pade1):
    t = np.linspace(0, 1, 2 * n + 1)
    res_t = np.zeros(len(t), dtype=complex)
    fit_t(scale * t, res_t, expn_pade1, etal_pade1)
    res_t = np.imag(res_t)
    print(res_t[-100:])
    n_sample = (len(t) + 1) // 2
    print(len(t), n_sample)
    h1 = res_t.copy()
    H = np.zeros((n_sample, n_sample))
    for i in range(n_sample):
        H[i, :] = h1[i:n_sample + i]

    vs, Q1 = LA.eigh(H)
    del H, t, res_t
    gc.collect()
    Q1 = Q1 * np.exp(-1j * np.angle(vs) / 2.0)
    vs = np.abs(vs)
    sort_array = np.argsort(vs)[::-1]
    vs = vs[sort_array]
    Q1 = (Q1[:, sort_array])
    del sort_array
    gc.collect()
    Q1 = (Q1[:, :20]).copy()
    gc.collect()
    return Q1, h1


def get_etal_expn(n, Q, n_gamma, h, etal_real, etar_real, expn_real):
    print("len of gamma", n_gamma)
    gamma = np.roots(Q[:, n_gamma][::-1])
    gamma = gamma[np.argsort(np.abs(gamma))[:n_gamma]]
    gamma = np.real(gamma)
    t_new = 2 * n * np.log(np.abs(gamma))

    gamma_m = np.zeros((n * 2 + 1, n_gamma))
    for i in range(n_gamma):
        for j in range(n * 2 + 1):
            gamma_m[j, i] = gamma[i]**j
    omega_new = np.dot(LA.inv(np.dot(np.transpose(gamma_m), gamma_m)),
                    np.dot(np.transpose(gamma_m), np.transpose(h)))

    etal = 1.j * omega_new
    etar = np.conjugate(1.j * omega_new)
    etaa = np.abs(omega_new)
    expn = -t_new / scale

    etal = np.append(etal, etal_real)
    etar = np.append(etar, etar_real)
    expn = np.append(expn, expn_real)
    return expn, etal, etar



def main():
    etal_pade1, etar_pade1, etaa_pade1, expn_pade1 = decompose_spe(
        phixx_sp, sp_para_dict, para_dict1, condition_dict, npsd, -1)
    etal_pade2, etar_pade2, etaa_pade2, expn_pade2 = decompose_spe(
        phiyy_sp, sp_para_dict, para_dict2, condition_dict, npsd, 1)

    etal_real1, etar_real1, etaa_real1, expn_real1 = decompose_spe_real(
        phixx_sp, sp_para_dict, para_dict1, condition_dict, 1)

    etal_real2, etar_real2, etaa_real2, expn_real2 = decompose_spe_real(
        phixx_sp, sp_para_dict, para_dict1, condition_dict, -1)

    n_spe = 1
    len_ = 10000
    spe_wid = 2000
    w = np.linspace(-spe_wid, spe_wid, len_)
    res_w = np.zeros(len(w), dtype=complex)

    phixx = lams1 * gams1**2 / (((w - mu_x) - omgs1)**2 + gams1**2)
    phixy = 0 * w
    phiyy = lams2 * gams2**2 / (((w - mu_x) - omgs2)**2 + gams2**2)

    fit_J(w, res_w, expn_pade1, etal_pade1, 1)
    plt.plot(w, (phixx / (1 + np.exp(beta * (w - mu_x)))) - res_w.real,
            'r',
            label='phixx')
    res_w = np.zeros(len(w), dtype=complex)
    fit_J(w, res_w, expn_pade2, etal_pade2, -1)
    plt.plot(w, (phiyy / (1 + np.exp(-beta * (w - mu_x)))) - res_w.real,
            'b',
            label='phixx')
    plt.legend(loc='best')
    plt.savefig('spe.pdf')
    plt.clf()
    del res_w, w, phixx, phixy, phiyy
    gc.collect()

    Q1, h1 = get_Q_h(n, expn_pade1, etal_pade1)
    Q2, h2 = get_Q_h(n, expn_pade2, etal_pade2)

    for npfs in npfs_l:
        expn1, etal1, etar1 = get_etal_expn(n, Q1, npfs, h1, etal_real1, etar_real1, expn_real1)
        expn2, etal2, etar2 = get_etal_expn(n,
        Q2, npfs, h2, etal_real2, etar_real2, expn_real2)

        expn_pade1 += -1.j * mu_x
        expn_pade2 += 1.j * mu_x
        expn1 += -1.j * mu_x
        expn2 += 1.j * mu_x

        np.savetxt("expn1_{}".format(npfs), expn1)
        np.savetxt("etal1_{}".format(npfs), etal1)
        np.savetxt("etar1_{}".format(npfs), etar1)
        np.savetxt("expn2_{}".format(npfs), expn2)
        np.savetxt("etal2_{}".format(npfs), etal2)
        np.savetxt("etar2_{}".format(npfs), etar2)

        len_ = 100000
        spe_wid = 200

        w = np.linspace(-spe_wid, spe_wid, len_)

        phixx = lams1 * gams1**2 / (((w - mu_x) - omgs1)**2 + gams1**2)
        phixy = 0 * w
        phiyy = lams2 * gams2**2 / (((w - mu_x) - omgs2)**2 + gams2**2)

        res_J1 = np.zeros(len(w), dtype=complex)
        res_J1_pade = np.zeros(len(w), dtype=complex)
        fit_J(w, res_J1, expn1, etal1, 1)
        fit_J(w, res_J1_pade, expn_pade1, etal_pade1, 1)

        res_J2 = np.zeros(len(w), dtype=complex)
        res_J2_pade = np.zeros(len(w), dtype=complex)
        fit_J(w, res_J2, expn2, etal2, -1)
        fit_J(w, res_J2_pade, expn_pade2, etal_pade2, -1)

        plt.plot(w, (phixx / (1 + np.exp(beta * (w - mu_x)))) - res_J1.real,
                'r',
                label='phixx')
        plt.plot(w, (phixx / (1 + np.exp(beta * (w - mu_x)))) - res_J1_pade.real,
                'b--',
                label='phixx')
        plt.plot(w, (phiyy / (1 + np.exp(-beta * (w - mu_x)))) - res_J2.real,
                'k',
                label='phixx')
        plt.plot(w, (phiyy / (1 + np.exp(-beta * (w - mu_x)))) - res_J2_pade.real,
                'g--',
                label='phixx')
        plt.savefig("error{}.pdf".format(npfs))
        plt.clf()

        len_ = 100000
        spe_wid = 1

        w = np.linspace(-spe_wid, spe_wid, len_)

        phixx = lams1 * gams1**2 / (((w - mu_x) - omgs1)**2 + gams1**2)
        phixy = 0 * w
        phiyy = lams2 * gams2**2 / (((w - mu_x) - omgs2)**2 + gams2**2)

        res_J1 = np.zeros(len(w), dtype=complex)
        res_J1_pade = np.zeros(len(w), dtype=complex)
        fit_J(w, res_J1, expn1, etal1, 1)
        fit_J(w, res_J1_pade, expn_pade1, etal_pade1, 1)

        res_J2 = np.zeros(len(w), dtype=complex)
        res_J2_pade = np.zeros(len(w), dtype=complex)
        fit_J(w, res_J2, expn2, etal2, -1)
        fit_J(w, res_J2_pade, expn_pade2, etal_pade2, -1)

        plt.plot(w, (phixx / (1 + np.exp(beta * (w - mu_x)))) - res_J1.real,
                'r',
                label='phixx')
        plt.plot(w, (phixx / (1 + np.exp(beta * (w - mu_x)))) - res_J1_pade.real,
                'b--',
                label='phixx')

        plt.plot(w, (phiyy / (1 + np.exp(-beta * (w - mu_x)))) - res_J2.real,
                'k',
                label='phixx')
        plt.plot(w, (phiyy / (1 + np.exp(-beta * (w - mu_x)))) - res_J2_pade.real,
                'g--',
                label='phixx')
        plt.savefig("error_zoom_in{}.pdf".format(npfs))

w_sp, lams1_sp, gams1_sp, omgs1_sp, gams2_sp, lams2_sp, omgs2_sp, beta_sp = sp.symbols(
    r"\omega , \lambda_1, \gamma_1, \Omega_1, \lambda_2, \gamma_2, \Omega_2, \beta",
    real=True)
phixx_sp = lams1_sp * gams1_sp**2 / ((w_sp - omgs1_sp)**2 + gams1_sp**2)
phixy_sp = 0
phiyy_sp = lams2_sp * gams2_sp**2 / ((w_sp - omgs2_sp)**2 + gams2_sp**2)

lmax = 10
dt = 0.01
n = 20000
scale = 3600
npsd = 1024
npfs_l = [9,11]
nmod = 4
lams1 = 1
lams2 = 1
gams1 = 1
gams2 = 1
omgs1 = 0
omgs2 = 0
temp = 0.001
beta = 1 / temp
mu_x = 0

sp_para_dict = {
    lams1_sp: lams1,
    lams2_sp: lams2,
    gams1_sp: gams1,
    gams2_sp: gams2,
    omgs1_sp: omgs1,
    omgs2_sp: omgs2
}

condition_dict = {}
para_dict1 = {'beta': beta}
para_dict2 = {'beta': beta}

main()
