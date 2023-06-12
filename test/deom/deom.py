import numpy as np
import sympy as sp
import os
from scipy.linalg import expm as matrix_exp
from deom.aux import fastread, PSD, fBose, direct_product_2d
import subprocess


def benchmark(file_str, magic_str, dir_str='bose', if_np=True, if_save=False):
    cmd = r'mv {} {}-{}'.format(file_str, file_str, magic_str)
    result = subprocess.call(cmd, shell=True)
    if os.path.exists('result/{}/{}-{}.npy'.format(dir_str, file_str, magic_str)):
        data1 = np.load(
            'result/{}/{}-{}.npy'.format(dir_str, file_str, magic_str))
    else:
        data1 = fastread_np(
            'result/{}/{}-{}'.format(dir_str, file_str, magic_str))
    data2 = fastread_np('./{}-{}'.format(file_str, magic_str))
    if if_save:
        np.save('result/{}/{}-{}'.format(dir_str, file_str, magic_str), data2)
    result = (np.sum(np.abs(data1 - data2)))
    if float(result) > 1e-6:
        print(result, 'FAILED')
    else:
        print(result, 'PASSED')


def thermal_equilibrium(beta, hams):
    return matrix_exp(- beta * hams) / np.trace(matrix_exp(- beta * hams))


def direct_product(a, *args):
    c = a.copy()
    for arg in args:
        if not isinstance(arg, np.ndarray):
            raise TypeError('Input must be numpy.ndarray')
        c = direct_product_2d(c, arg)
    return c


def decompose_spe(spe, w_sp, sp_para_dict, para_dict, condition_dict, npsd, pade=1):
    if (sp.cancel(
            spe.subs(condition_dict)).as_real_imag()[1] == 0):
        imag_part = sp.cancel(
            spe.subs(condition_dict)).as_real_imag()[0]
    else:
        imag_part = sp.cancel(
            spe.subs(condition_dict)).as_real_imag()[1]
    numer, denom = sp.cancel(sp.factor(imag_part)).as_numer_denom()
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
    expn_val_n_cc = expn[expn_imag_sort[expn_imag == 0]]

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


def fastread_np(str):
    return fastread(str).to_numpy()


def sum_exptontial_spectrum(w, res, expn, etal, sigma):
    for i in range(len(etal)):
        res += etal[i] / (expn[i] + sigma * 1.j * w)


def convert(o):
    if isinstance(o, np.int64):
        return int(o)
    elif isinstance(o, np.int32):
        return int(o)
    raise TypeError


def sort_symmetry(etal, expn, if_sqrt=True):
    expn_imag_sort = np.argsort(np.abs(np.imag(expn)))[::-1]
    expn_imag = np.sort(np.abs(np.imag(expn)))[::-1]
    expn = expn[expn_imag_sort]
    etal = etal[expn_imag_sort]
    etar = etal[expn_imag_sort]
    expn_val_cc = np.where(expn[expn_imag > 1e-10])[0]
    etaa = np.zeros(len(etal), dtype=float)
    for ii in range(0, len(expn_val_cc), 2):
        even_i = ii
        odd_i = ii + 1
        etar[even_i] = np.conj(etal[odd_i])
        etar[odd_i] = np.conj(etal[even_i])
        etaa[even_i] = np.abs(etal[even_i])
        etaa[odd_i] = np.abs(etal[odd_i])
    for ii in range(len(expn_val_cc), len(expn)):
        even_i = ii
        etar[even_i] = np.conj(etal[even_i])
        etaa[even_i] = np.abs(etal[even_i])
    if (if_sqrt):
        etaa = np.sqrt(etaa)
    return expn, etal, etar, etaa


def complex_2_json(list_input):
    if (type(list_input) == np.ndarray):
        return {
            "if_initial": True,
            "real": list(np.real(list_input.flatten())),
            "imag": list(np.imag(list_input.flatten()))
        }
    else:
        return {
            "real": np.real(list_input),
            "imag": np.imag(list_input)
        }


def init_qmd(json_init, qmd1a, qmd1c, mode, nsys, etaa, etal, etar):
    qmdta_l = np.zeros((len(mode), nsys, nsys), dtype=complex)
    qmdta_r = np.zeros((len(mode), nsys, nsys), dtype=complex)
    qmdtc_l = np.zeros((len(mode), nsys, nsys), dtype=complex)
    qmdtc_r = np.zeros((len(mode), nsys, nsys), dtype=complex)
    for i in range(len(mode)):
        i_mod = mode[i]
        qmdta_l[i, :, :] = qmd1a[i_mod, :, :] * np.sqrt(etaa[i])
        qmdta_r[i, :, :] = qmd1a[i_mod, :, :] * np.sqrt(etaa[i])
        qmdtc_l[i, :, :] = qmd1c[i_mod, :, :] * etal[i] / np.sqrt(etaa[i])
        qmdtc_r[i, :, :] = qmd1c[i_mod, :, :] * etar[i] / np.sqrt(etaa[i])
    json_init["qmdta_l"] = complex_2_json(qmdta_l)
    json_init["qmdta_r"] = complex_2_json(qmdta_r)
    json_init["qmdtc_l"] = complex_2_json(qmdtc_l)
    json_init["qmdtc_r"] = complex_2_json(qmdtc_r)


# Do some normalize thing, you can find more details in the pdf file.
def init_qmd_quad(json_init, qmd2a, qmd2b, qmd2c, mode, nsys, nind, nmod, etaa, etal, etar):
    qmdt2a_l = np.zeros((nind*nind, nsys, nsys), dtype=complex)
    qmdt2a_r = np.zeros((nind*nind, nsys, nsys), dtype=complex)
    qmdt2b_l = np.zeros((nind*nind, nsys, nsys), dtype=complex)
    qmdt2b_r = np.zeros((nind*nind, nsys, nsys), dtype=complex)
    qmdt2c_l = np.zeros((nind*nind, nsys, nsys), dtype=complex)
    qmdt2c_r = np.zeros((nind*nind, nsys, nsys), dtype=complex)
    for i in range(len(mode)):
        for j in range(len(mode)):
            i_mod = mode[i]
            j_mod = mode[j]
            index_mat = i * nind + j
            qmdt2a_l[index_mat, :, :] = qmd2a[i_mod, j_mod, :, :] * \
                np.sqrt(etaa[i]) * np.sqrt(etaa[j])
            qmdt2a_r[index_mat, :, :] = qmd2a[i_mod, j_mod, :, :] * \
                np.sqrt(etaa[i]) * np.sqrt(etaa[j])
            qmdt2b_l[index_mat, :, :] = qmd2b[i_mod, j_mod, :, :] * \
                etal[i] / np.sqrt(etaa[i]) * np.sqrt(etaa[j])
            qmdt2b_r[index_mat, :, :] = qmd2b[i_mod, j_mod, :, :] * \
                etar[i] / np.sqrt(etaa[i]) * np.sqrt(etaa[j])
            qmdt2c_l[index_mat, :, :] = qmd2c[i_mod, j_mod, :, :] * \
                etal[i] * etal[j] / np.sqrt(etaa[i]) / np.sqrt(etaa[j])
            qmdt2c_r[index_mat, :, :] = qmd2c[i_mod, j_mod, :, :] * \
                etar[i] * etar[j] / np.sqrt(etaa[i]) / np.sqrt(etaa[j])
    json_init["qmdt2a_l"] = complex_2_json(qmdt2a_l)
    json_init["qmdt2a_r"] = complex_2_json(qmdt2a_r)
    json_init["qmdt2b_l"] = complex_2_json(qmdt2b_l)
    json_init["qmdt2b_r"] = complex_2_json(qmdt2b_r)
    json_init["qmdt2c_l"] = complex_2_json(qmdt2c_l)
    json_init["qmdt2c_r"] = complex_2_json(qmdt2c_r)
