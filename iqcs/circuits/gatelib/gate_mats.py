import numpy as np
from typing  import Union,  List
from scipy.linalg import sqrtm

def constarr(mat):
    mat.setflags(write=False)
    return mat

ID_MAT = constarr(np.asarray([[1., 0.], 
                              [0., 1.]], dtype=complex))

H_MAT = constarr(np.asarray([[1. / np.sqrt(2), 1. / np.sqrt(2)],
                             [1. / np.sqrt(2), -1. / np.sqrt(2)]], dtype=complex))

X_MAT = constarr(np.asarray([[0., 1.],
                             [1., 0.]], dtype=complex))

Y_MAT = constarr(np.asarray([[0.0, -1.0j], 
                             [1.0j, 0.0]], dtype=complex))

Z_MAT = constarr(np.asarray([[1.0, 0.0], 
                             [0.0, -1.0]], dtype=complex))

S_MAT = constarr(np.asarray([[1.0, 0.0], 
                           [0.0, 1.0j]], dtype=complex))


SX_MAT = constarr(sqrtm(X_MAT))

SY_MAT = constarr(sqrtm(Y_MAT))

T_MAT = constarr(np.asarray([[1.0, 0.0], 
                             [0.0, np.exp(1.0j * np.pi / 4)]], dtype=complex))

W_MAT = constarr((X_MAT + Y_MAT) / np.sqrt(2))

SW_MAT = constarr(sqrtm(W_MAT))

H_MAT = constarr((X_MAT + Z_MAT) / np.sqrt(2))

CX_MAT = constarr(np.asarray(
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
        [0.0, 0.0, 1.0, 0.0],
    ],
    dtype=complex,
))

SWAP_MAT =  constarr(np.asarray(
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ],
    dtype=complex,
))

ISWAP_MAT = constarr(np.asarray(
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 1.0j, 0.0],
        [0.0, 1.0j, 0.0, 0.0],
        [0.0, 0.0, 0.0, 1.0],
    ],
    dtype=complex,
))

CY_MAT = constarr(np.asarray(
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, -1.0j],
        [0.0, 0.0, 1.0j, 0.0],
    ],
    dtype=complex,
))

CZ_MAT = constarr(np.asarray(
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, -1.0],
    ],
    dtype=complex,
))

CS_MAT = constarr(np.asarray(
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, 1.0j],
    ],
    dtype=complex,
))

CT_MAT = constarr(np.asarray(
    [
        [1.0, 0.0, 0.0, 0.0],
        [0.0, 1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0, 0.0],
        [0.0, 0.0, 0.0, np.exp(1.0j * np.pi / 4)],
    ],
    dtype=complex,
))

def u2matrix(paras):
    "OpenQASM 3.0 specification"
    _phi, _lambda = paras
    return np.array([[1., np.exp(-1.j * _lambda)],
                     [np.exp(1.j * _phi), np.exp((_phi + _lambda) * 1.j)]], dtype=complex)


def u3matrix(paras):
    "OpenQASM 3.0 specification"
    _theta, _phi, _lambda = paras
    return np.array([[np.cos(0.5 * _theta), -np.exp(_lambda * 1.j) * np.sin(0.5 * _theta)],
                     [np.exp(_phi * 1.j) * np.sin(0.5 * _theta),
                      np.exp((_phi + _lambda) * 1.j) * np.cos(0.5 * _theta)]], dtype=complex)


def rxmatrix(para):
    theta, = para
    return np.array([[np.cos(0.5 * theta), -1.j * np.sin(0.5 * theta)],
                     [-1.j * np.sin(0.5 * theta), np.cos(0.5 * theta)]], dtype=complex)


def rymatrix(para):
    theta, = para
    return np.array([[np.cos(0.5 * theta), - np.sin(0.5 * theta)],
                     [np.sin(0.5 * theta), np.cos(0.5 * theta)]], dtype=complex)


def rzmatrix(para):
    theta, = para
    return np.array([[np.exp(-0.5j * theta), 0.],
                     [0., np.exp(0.5j * theta)]], dtype=complex)

def pmatrix(para):
    labda, = para
    return np.array([[1, 0], 
                     [0, np.exp(1j*labda)]] ,dtype=complex)

def rxxmatrix(para):
    """Unitary evolution of XX interaction"""
    theta, = para
    return np.array([[np.cos(theta/2), 0, 0, -1j*np.sin(theta/2)],
                     [0, np.cos(theta/2), -1j*np.sin(theta/2), 0],
                     [0, -1j*np.sin(theta/2), np.cos(theta/2), 0],
                     [-1j*np.sin(theta/2), 0, 0, np.cos(theta/2)]
                    ])

def ryymatrix(para):
    """ Unitary evolution of YY interaction"""
    theta, = para
    return np.array([[np.cos(theta/2), 0, 0, 1j*np.sin(theta/2)],
                     [0, np.cos(theta/2), -1j*np.sin(theta/2), 0],
                     [0, -1j*np.sin(theta/2), np.cos(theta/2), 0],
                     [1j*np.sin(theta/2), 0, 0, np.cos(theta/2)]
                    ])

def rzzmatrix(para):
    theta, = para
    return np.array([[np.exp(-1j*theta/2), 0, 0, 0],
                     [0, np.exp(1j*theta/2), 0, 0],
                     [0, 0, np.exp(1j*theta/2), 0],
                     [0, 0, 0, np.exp(-1j*theta/2)]
                    ])