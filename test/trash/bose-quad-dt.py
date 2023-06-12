import math
import os
import time
import numpy as np
import sympy as sp
import subprocess
import itertools
import json
import matplotlib.pyplot as plt
from cvxopt import solvers, matrix, spmatrix, mul
from numpy import linalg as LA
from scipy import sparse
import pandas as pd
from termcolor import colored


def fastread(str):
    try:
        data = pd.read_csv(str, header=None, sep='\s+', dtype=np.float64)
    except:
        data = pd.read_csv(str,
                           header=None,
                           sep='\s+',
                           skipfooter=1,
                           dtype=np.float64)
    return data


def fastread_np(str):
    return fastread(str).to_numpy()


def convert(o):
    if isinstance(o, np.int64):
        return int(o)
    raise TypeError


def fBose(x, pole, resi):
    return 1 / x + 0.5 + sum(2.0 * resi[i] * x / (x**2 + pole[i]**2)
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


def PSD(N, BoseFermi=1, pade=1):
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


def arma_print(ndarray):

    shape = ndarray.shape
    dimen = len(shape)

    if dimen == 1:

        if issubclass(type(ndarray[0]), np.int_):
            print('ARMA_MAT_TXT_IS004\n%d %d' % (shape[0], 1))
            for row in ndarray:
                print('%d' % row)
        elif issubclass(type(ndarray[0]), float):
            print('ARMA_MAT_TXT_FN008\n%d %d' % (shape[0], 1))
            for row in ndarray:
                print('%.8e' % row)
        elif issubclass(type(ndarray[0]), complex):
            print('ARMA_MAT_TXT_FC016\n%d %d' % (shape[0], 1))
            for row in ndarray:
                print('(%.8e,%-.8e)' % (row.real, row.imag))

    elif dimen == 2:

        if issubclass(type(ndarray[0, 0]), np.int_):
            print('ARMA_MAT_TXT_IS004\n%d %d' % (shape[0], shape[1]))
            for row in ndarray:
                print(' '.join('%d' % x for x in row))
        elif issubclass(type(ndarray[0, 0]), float):
            print('ARMA_MAT_TXT_FN008\n%d %d' % (shape[0], shape[1]))
            for row in ndarray:
                print(' '.join('%.8e' % x for x in row))
        elif issubclass(type(ndarray[0, 0]), complex):
            print('ARMA_MAT_TXT_FC016\n%d %d' % (shape[0], shape[1]))
            for row in ndarray:
                print(' '.join('(%.8e,%-.8e)' % (x.real, x.imag) for x in row))

    elif dimen == 3:

        if issubclass(type(ndarray[0, 0, 0]), np.int_):
            print('ARMA_CUB_TXT_IS004\n%d %d %d' %
                  (shape[1], shape[2], shape[0]))
            for slc in ndarray:
                for row in slc:
                    print(' '.join('%d' % x for x in row))
        elif issubclass(type(ndarray[0, 0, 0]), float):
            print('ARMA_CUB_TXT_FN008\n%d %d %d' %
                  (shape[1], shape[2], shape[0]))
            for slc in ndarray:
                for row in slc:
                    print(' '.join('%-.8e' % x for x in row))
        elif issubclass(type(ndarray[0, 0, 0]), complex):
            print('ARMA_CUB_TXT_FC016\n%d %d %d' %
                  (shape[1], shape[2], shape[0]))
            for slc in ndarray:
                for row in slc:
                    print(' '.join('(%.8e,%-.8e)' % (x.real, x.imag)
                                   for x in row))


def arma_write(ndarray, filename):

    shape = ndarray.shape
    dimen = len(shape)

    with open(filename, 'w') as f:
        if dimen == 1:
            if issubclass(type(ndarray[0]), np.int_):
                print('ARMA_MAT_TXT_IS004\n%d %d' % (shape[0], 1), file=f)
                for row in ndarray:
                    print('%d' % row, file=f)
            elif issubclass(type(ndarray[0]), float):
                print('ARMA_MAT_TXT_FN008\n%d %d' % (shape[0], 1), file=f)
                for row in ndarray:
                    print('%.8e' % row, file=f)
            elif issubclass(type(ndarray[0]), complex):
                print('ARMA_MAT_TXT_FC016\n%d %d' % (shape[0], 1), file=f)
                for row in ndarray:
                    print('(%.8e,%-.8e)' % (row.real, row.imag), file=f)

        elif dimen == 2:

            if issubclass(type(ndarray[0, 0]), np.int_):
                print('ARMA_MAT_TXT_IS004\n%d %d' % (shape[0], shape[1]),
                      file=f)
                for row in ndarray:
                    print(' '.join('%d' % x for x in row), file=f)
            elif issubclass(type(ndarray[0, 0]), float):
                print('ARMA_MAT_TXT_FN008\n%d %d' % (shape[0], shape[1]),
                      file=f)
                for row in ndarray:
                    print(' '.join('%.8e' % x for x in row), file=f)
            elif issubclass(type(ndarray[0, 0]), complex):
                print('ARMA_MAT_TXT_FC016\n%d %d' % (shape[0], shape[1]),
                      file=f)
                for row in ndarray:
                    print(' '.join('(%.8e,%-.8e)' % (x.real, x.imag)
                                   for x in row),
                          file=f)

        elif dimen == 3:

            if issubclass(type(ndarray[0, 0, 0]), np.int_):
                print('ARMA_CUB_TXT_IS004\n%d %d %d' %
                      (shape[1], shape[2], shape[0]),
                      file=f)
                for slc in ndarray:
                    for row in slc:
                        print(' '.join('%d' % x for x in row))
            elif issubclass(type(ndarray[0, 0, 0]), float):
                print('ARMA_CUB_TXT_FN008\n%d %d %d' %
                      (shape[1], shape[2], shape[0]),
                      file=f)
                for slc in ndarray:
                    for row in slc:
                        print(' '.join('%-.8e' % x for x in row), file=f)
            elif issubclass(type(ndarray[0, 0, 0]), complex):
                print('ARMA_CUB_TXT_FC016\n%d %d %d' %
                      (shape[1], shape[2], shape[0]),
                      file=f)
                for slc in ndarray:
                    for row in slc:
                        print(' '.join('(%.8e,%-.8e)' % (x.real, x.imag)
                                       for x in row),
                              file=f)


# in this script, we can decompose any given spectrum, but the sympy format is must been given
# do u like haskell?
# sympy[spe(def by sympy)], dict[sp_para_dict], dict[para_dict], dict[npsd],
# dict[pade] >> np.array[etal], np.array[etar],np.array[etaa], np.array[expn]


def decompose_spe(spe, sp_para_dict, para_dict, condition_dict, npsd, pade=1):
    numer, denom = sp.cancel(
        sp.factor(sp.cancel(
            spe.subs(condition_dict)).as_real_imag()[1])).as_numer_denom()
    numer_get_para = (sp.factor(numer)).subs(sp_para_dict)
    denom_get_para = (sp.factor(denom)).subs(sp_para_dict)

    poles = sp.nroots(denom_get_para)
    float(sp.re(poles[0]))

    expn = []
    poles_allplane = np.array([])
    for i in poles:
        i = complex(i)
        if i.imag < 0:
            expn.append(i * 1.J)
        poles_allplane = np.append(poles_allplane, i)

    etal = []
    etar = []
    etaa = []

    expn = np.array(expn)

    expn_imag_sort = np.argsort(np.abs(np.imag(expn)))[::-1]
    expn_imag = np.sort(np.abs(np.imag(expn)))[::-1]

    expn_val_cc = expn[expn_imag_sort[expn_imag != 0]]
    # expn_arg_cc = expn_imag_sort[expn_imag != 0]
    expn_val_n_cc = expn[expn_imag_sort[expn_imag == 0]]
    # expn_arg_n_cc = expn_imag_sort[expn_imag == 0]

    expn = list(expn[expn_imag_sort])
    pole, resi = PSD(npsd, 1, pade)
    beta = para_dict['beta']
    temp = 1 / beta

    for ii in range(0, len(expn_val_cc), 2):
        etal.append(
            complex(
                sp.N((-2.j * numer_get_para /
                      np.multiply.reduce(w_sp - poles_allplane[np.abs(
                          poles_allplane + 1.J * expn_val_cc[ii]) > 1e-14])
                      ).subs({w_sp: -1.j * expn_val_cc[ii]}) *
                     fBose(-1.J * expn_val_cc[ii] / temp, pole, resi))))

        etal.append(
            complex(
                sp.N((-2.j * numer_get_para /
                      np.multiply.reduce(w_sp - poles_allplane[np.abs(
                          poles_allplane + 1.J * expn_val_cc[ii + 1]) > 1e-14])
                      ).subs({w_sp: -1.j * expn_val_cc[ii + 1]}) *
                     fBose(-1.J * expn_val_cc[ii + 1] / temp, pole, resi))))

        etar.append(np.conj(etal[-1]))
        etar.append(np.conj(etal[-2]))
        etaa.append(np.sqrt(np.abs(etal[-2]) * np.abs(etar[-2])))
        etaa.append(np.sqrt(np.abs(etal[-1]) * np.abs(etar[-1])))

    for ii in range(len(expn_val_n_cc)):
        etal.append(
            complex(
                sp.N((-2.j * numer_get_para /
                      np.multiply.reduce(w_sp - poles_allplane[np.abs(
                          poles_allplane + 1.J * expn_val_n_cc[ii]) > 1e-14])
                      ).subs({w_sp: -1.j * expn_val_n_cc[ii]}) *
                     fBose(-1.J * expn_val_n_cc[ii] / temp, pole, resi))))
        etar.append(np.conj(etal[-1]))
        etaa.append(np.sqrt(np.abs(etal[-1]) * np.abs(etar[-1])))

    f = numer_get_para / np.multiply.reduce(w_sp - poles_allplane)
    f = sp.lambdify(w_sp, f)

    for inma in range(len(pole)):
        zomg = -1.J * pole[inma] * temp
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


def fit_J(w, expn, etal):
    res = np.zeros_like(w, dtype=complex)
    for i in range(len(etal)):
        res += etal[i] / (expn[i] - 1.j * w)
    return res


def INDEX3(i, j, k, mum):
    return mum * mum * i + mum * j + k


def fit_t(t, res, expn, etal):
    for i in range(len(etal)):
        res += etal[i] * np.exp(-expn[i] * t)
    return res


def numpy_to_cvxopt_matrix(A):
    if A is None:
        return A
    elif isinstance(A, np.ndarray):
        if A.ndim == 1:
            return matrix(A, (A.shape[0], 1), 'd')
        else:
            return matrix(A, A.shape, 'd')
    else:
        return A


# for vib_key, S_0, npsd, lmax in itertools.product([0, 0.5, 1.0], [0.25, 0.5], [1, 2], [40, 50, 60]):
# para:2018 jcp 114103
for npsd, lmax, (theta, eta), init, ferr, dt in itertools.product(
    [1],  # npsd
    [100],  # lmax
    [(3 / 4, 0.25)],  # theta, eta
    [1],  # init, 1 for absorption, 0 for emission
    [1e-10],  # ferr
    [0.01],  # dt
):
    nmod = 1
    weg = 1
    temp = 1
    OMG = 1
    beta = 1 / temp
    omgb = 1
    ti_gam = 15
    ti_eta = 10
    sdip_val, bdip0_val, bdip1_val, bdip2_val = 1, 1, 1, 1
    if (init == 1):
        omgb_BO = 1
    if (init == 0):
        omgb_BO = theta

    if (init == 1):
        alp0 = eta * theta * theta
        alp1 = -np.sqrt(2 * eta * omgb) * theta * theta
        alp2 = omgb * (theta * theta - 1) / 2
    if (init == 0):
        alp0 = eta
        alp1 = np.sqrt(2 * eta * omgb / theta)
        alp2 = -omgb * (theta * theta - 1) / 2 / theta

    print(alp0, alp1, alp2)

    w_sp, omgs_sp, omgs_p_sp, eta_sp, gamma_sp, beta_sp = sp.symbols(
        r"\omega , \omega_{b}, \omega'_{b}, \eta, \gamma, \beta", real=True)
    phixx_sp = omgs_p_sp / (
        (omgs_p_sp * omgs_p_sp - w_sp * w_sp) + w_sp * eta_sp * omgs_sp /
        (w_sp + sp.I * gamma_sp))
    spe_vib_sp = phixx_sp
    sp_para_dict = {
        omgs_p_sp: omgb_BO,
        omgs_sp: omgb,
        eta_sp: ti_eta,
        gamma_sp: ti_gam
    }
    condition_dict = {}
    para_dict = {'beta': beta}
    etal_pade, etar_pade, etaa_pade, expn_pade = decompose_spe(
        spe_vib_sp, sp_para_dict, para_dict, condition_dict, 200)
    n = 2500
    scale = 50
    t = np.linspace(0, 1, 2 * n + 1)
    res_t = np.zeros(len(t), dtype=complex)
    fit_t(scale * t, res_t, expn_pade, etal_pade)
    n_sample = (len(t) + 1) // 2
    h = res_t
    H = np.zeros((n_sample, n_sample))
    for i in range(n_sample):
        H[i, :] = np.imag(h[i:n_sample + i])
    sing_vs, Q = LA.eig(H)
    phase_mat = np.diag(
        [np.exp(-1j * np.angle(sing_v) / 2.0) for sing_v in sing_vs])
    vs = np.array([np.abs(sing_v) for sing_v in sing_vs])
    Qp = np.dot(Q, phase_mat)
    sort_array = np.argsort(vs)[::-1]
    vs = vs[sort_array]
    Qp = (Qp[:, sort_array])
    print(vs[:20])
    for n_gamma in [2]:
        print("len of gamma", n_gamma)
        for i in [n_gamma]:
            print(i)
            gamma = np.roots(Qp[:, i][::-1])
        gamma_new = gamma[np.argsort(np.abs(gamma))[:n_gamma]]
        t_new = 2 * n * np.log(gamma_new)
    n_col = n_sample * 2 - 1
    n_row = n_gamma
    gamma_m = np.zeros((2 * n_col, 2 * n_row), dtype=float)
    for i in range(n_row):
        for j in range(n_col):
            gamma_m[j, i] = np.real(gamma_new[i]**j)
            gamma_m[n_col + j, n_row + i] = np.real(gamma_new[i]**j)
            gamma_m[j, n_row + i] = -np.imag(gamma_new[i]**j)
            gamma_m[n_col + j, i] = np.imag(gamma_new[i]**j)
    h_m = np.append(np.real(h), np.imag(h))
    freq_d = np.append(
        np.append(np.linspace(-10000, 10, n_col // 2),
                  np.linspace(-10, 10, n_col + 1)),
        np.linspace(10, 10000, n_col // 2))

    freq_m = np.zeros((2 * n_col, 2 * n_row), dtype=float)
    expn = -t_new / scale
    for i in range(n_row):
        for j in range(2 * n_col):
            freq_m[j,
                   i] = np.real(expn[i]) / (np.real(expn[i])**2 +
                                            (np.imag(expn[i]) - freq_d[j])**2)
            freq_m[j, n_row + i] = (np.imag(expn[i]) - freq_d[j]) / \
                (np.real(expn[i])**2 + (np.imag(expn[i]) - freq_d[j])**2)
    C = numpy_to_cvxopt_matrix(gamma_m)
    d = numpy_to_cvxopt_matrix(h_m)
    A = numpy_to_cvxopt_matrix(-freq_m)
    b = numpy_to_cvxopt_matrix(np.zeros(2 * n_col))
    Q = C.T * C
    q = -d.T * C
    opts = {
        'show_progress': False,
        'abstol': 1e-24,
        'reltol': 1e-24,
        'feastol': 1e-24
    }
    for k, v in opts.items():
        solvers.options[k] = v
    sol = solvers.qp(Q, q.T, A, b, None, None, None, None)
    omega_new_temp = np.array(sol['x']).reshape(2, n_gamma)
    omega_new = omega_new_temp[0, :] + 1.j * omega_new_temp[1, :]
    etal = omega_new
    expn = -t_new / scale
    expn_imag_sort = np.argsort(np.abs(np.imag(expn)))[::-1]
    expn_imag = np.sort(np.abs(np.imag(expn)))[::-1]
    expn = expn[expn_imag_sort]
    etal = etal[expn_imag_sort]
    etar = etal[expn_imag_sort]
    expn_val_cc = np.where(expn[expn_imag > 1e-14])[0]
    etaa = np.zeros(len(etal), dtype=float)
    for ii in range(0, len(expn_val_cc), 2):
        even_i = ii
        odd_i = ii + 1
        imag_mean = (np.imag(etal[even_i]) - np.imag(etal[odd_i])) / 2
        etal[even_i] = np.real(etal[even_i]) + 1.j * imag_mean
        etal[odd_i] = np.real(etal[odd_i]) - 1.j * imag_mean
        etar[even_i] = np.conj(etal[odd_i])
        etar[odd_i] = np.conj(etal[even_i])
        etaa[even_i] = np.abs(etal[even_i])
        etaa[odd_i] = np.abs(etal[odd_i])
    for ii in range(len(expn_val_cc), n_gamma):
        even_i = ii
        etar[even_i] = np.conj(etal[even_i])
        etaa[even_i] = np.abs(etal[even_i])

    len_ = 100000
    spe_wid = 20
    w = np.append(np.linspace(0, spe_wid, len_),
                  np.linspace(-spe_wid, 0, len_))
    phixx = np.imag(omgb_BO /
                    ((omgb_BO * omgb_BO - w * w) + w * ti_eta * omgb /
                     (w + 1.j * ti_gam)))
    spe_vib = phixx

    plt.plot(w[:len_], (spe_vib / (1 - np.exp(-beta * w)))[:len_] -
             fit_J(w, expn, etal).real[:len_], 'r')
    plt.plot(w[len_:], (spe_vib / (1 - np.exp(-beta * w)))[len_:] -
             fit_J(w, expn, etal).real[len_:], 'r')
    plt.savefig('spe_{}_{}.pdf'.format(str(lmax), str(npsd)))
    plt.clf()

    mode = np.zeros_like(expn, dtype=int)

    # syst \e> as 0 |g>as 1
    hams = np.zeros((2, 2), dtype=complex)
    hams[0, 0] = 1
    hams[0, 1] = 1
    hams[1, 0] = 1

    qmds = np.zeros((nmod, 2, 2), dtype=complex)
    if (init == 1):
        qmds[0, 0, 0] = 1
    if (init == 0):
        qmds[0, 1, 1] = 1

    # hidx
    trun = 0
    nmax = 200000

    sdip = np.zeros((nmod, 2, 2), dtype=float)
    sdip[0, 0, 1] = sdip_val * 1.0
    sdip[0, 1, 0] = sdip_val * 1.0

    pdip0 = np.zeros((nmod, 2, 2), dtype=float)
    pdip1 = np.zeros((nmod, 2, 2), dtype=float)
    pdip2 = np.zeros((nmod, 2, 2), dtype=float)

    bdip0 = bdip0_val * np.zeros(12 + 2 * npsd, dtype=float)
    bdip1 = bdip1_val * np.zeros(12 + 2 * npsd, dtype=float)
    bdip2 = bdip2_val * np.zeros(12 + 2 * npsd, dtype=float)

    sdipN = np.zeros((nmod, 2, 2), dtype=float)
    pdipN = np.zeros((nmod, 2, 2), dtype=float)
    bdipN = np.zeros(12 + 2 * npsd, dtype=float)

    # proprho
    jsonInit = {
        "nmax": nmax,
        "lmax": lmax,
        "alp0": alp0,
        "alp1": alp1,
        "alp2": alp2,
        "ferr": ferr,
        "nind": len(expn),
        "nmod": nmod,
        "filter": False,
        "equilibrium": {
            "sc2": False,
            "cgsb": False,
            "dt-method": True,
            "ti": 0,
            "tf": 10,
            "dt": dt,
            "der": False,
            "OMG": OMG,
            "backup": True,
        },
        "corr": {
            "ti": 0,
            "tf": 0,
            "dt": dt,
            "re_tree": 1,
            "filter": 1,
        },
        "omgb1": {
            "real": np.real(omgb_BO),
            "imag": np.imag(omgb_BO)
        },
        "expn": {
            "real": list(np.real(expn.flatten())),
            "imag": list(np.imag(expn.flatten()))
        },
        "ham1": {
            "real": list(np.real(hams.flatten())),
            "imag": list(np.imag(hams.flatten()))
        },
        "qmd1": {
            "real": list(np.real(qmds.flatten())),
            "imag": list(np.imag(qmds.flatten()))
        },
        "modLabel": {
            "real": list(np.real(mode.flatten())),
            "imag": list(np.imag(mode.flatten()))
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
                "real": list(np.real(sdip.flatten())),
                "imag": list(np.imag(sdip.flatten()))
            },
            "bdip0_cub": {
                "real": list(np.real(bdip0.flatten())),
                "imag": list(np.imag(bdip0.flatten()))
            },
            "pdip0_cub": {
                "real": list(np.real(pdip0.flatten())),
                "imag": list(np.imag(pdip0.flatten()))
            },
            "bdip1_cub": {
                "real": list(np.real(bdip1.flatten())),
                "imag": list(np.imag(bdip1.flatten()))
            },
            "pdip1_cub": {
                "real": list(np.real(pdip1.flatten())),
                "imag": list(np.imag(pdip1.flatten()))
            },
            "bdip2_cub": {
                "real": list(np.real(bdip2.flatten())),
                "imag": list(np.imag(bdip2.flatten()))
            },
            "pdip2_cub": {
                "real": list(np.real(pdip2.flatten())),
                "imag": list(np.imag(pdip2.flatten()))
            }
        },
        "dipole1": {
            "sdip_cub": {
                "real": list(np.real(sdip.flatten())),
                "imag": list(np.imag(sdip.flatten()))
            },
            "bdip0_cub": {
                "real": list(np.real(bdip0.flatten())),
                "imag": list(np.imag(bdip0.flatten()))
            },
            "pdip0_cub": {
                "real": list(np.real(pdip0.flatten())),
                "imag": list(np.imag(pdip0.flatten()))
            },
            "bdip1_cub": {
                "real": list(np.real(bdip1.flatten())),
                "imag": list(np.imag(bdip1.flatten()))
            },
            "pdip1_cub": {
                "real": list(np.real(pdip1.flatten())),
                "imag": list(np.imag(pdip1.flatten()))
            },
            "bdip2_cub": {
                "real": list(np.real(bdip2.flatten())),
                "imag": list(np.imag(bdip2.flatten()))
            },
            "pdip2_cub": {
                "real": list(np.real(pdip2.flatten())),
                "imag": list(np.imag(pdip2.flatten()))
            }
        },
    }

    with open('input.json', 'w') as f:
        json.dump(jsonInit, f, indent=4, default=convert)
    cmd = r'export OMP_NUM_THREADS={}'.format(8)
    cmd += '&&' + r'/usr/local/bin/jemalloc.sh  ../bose_quad_2.out'
    start_time = time.time()
    with open('out-equilibrium', "w") as outfile:
        result = subprocess.call(cmd, shell=True, stdout=outfile)
    np.savetxt('time-equilibrium', [time.time() - start_time])
    cmd = r'mv prop-rho-eq1.dat prop-rho-eq1.dat-equilibrium'
    os.system(cmd)
    data1 = fastread_np('./result/prop-rho-eq1.dat-equilibrium')
    data2 = fastread_np('./prop-rho-eq1.dat-equilibrium')
    result = (np.sum(np.abs(data1 - data2)))
    print(
        result,
        colored('PASSED', 'green') if float(result) < 1e-4 else colored(
            'FAILED', 'red'))
