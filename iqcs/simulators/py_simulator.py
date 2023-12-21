#default circuit simulator for state vector

from typing import Iterable, List, Union
from iqcs.circuits.quantum_circuit import QuantumCircuit
from ..circuits.gates import QuantumGate
from ..circuits.operations import  Delay, Barrier
import numpy as np
from sparse import COO, SparseArray
from scipy.sparse import kron, eye, coo_matrix
from ..utils.paulis import PauliOp

import copy

def global_op(gate: QuantumGate, global_qubits: List) -> coo_matrix:
    """Local operators to global operators"""
    num  = len(global_qubits)
    if len(gate.pos == 1):
        local_mat = coo_matrix(gate.matrix)     
        pos = global_qubits.index(gate.pos)
        local_mat = kron(kron(eye(2**pos), local_mat), eye(2**(num - pos-1)))
        return local_mat

    else:
        local_mat =coo_matrix(gate.matrix)
        pos = [global_qubits.index(p) for p in gate.pos]
        num_left = min(pos)
        num_right = num - max(pos) - 1
        num_center = max(pos) - min(pos) + 1
        center_mat = kron(local_mat, eye(2**(num_center - len(pos))))
        origin_order = sorted(pos)
        origin_order.extend([p for p in range(min(pos), max(pos)+1) if p not in pos])
        new_order = np.argsort(origin_order)
        center_mat = COO.from_scipy_sparse(center_mat)
        center_mat = permutebits(center_mat, new_order).to_scipy_sparse()
        center_mat = kron(kron(eye(2**num_left), center_mat), eye(2**num_right))
        return center_mat


def permutebits(mat: Union[SparseArray, np.ndarray], order : Iterable)\
    ->Union[SparseArray, np.ndarray]:
    """permute qubits for operators or states"""
    num = len(order)
    order = np.array(order)
    r = len(mat.shape)
    mat = np.reshape(mat, [2]*r*num)
    order = np.concatenate([order+len(order)*i for i  in range(r)]) 
    mat = np.transpose(mat, order)
    mat = np.reshape(mat, [2**num]*r)
    return mat

def ptrace(psi, ind_A: List, diag: bool=True) -> np.ndarray:
    """partial trace on a state vector, only for big endian form"""
    num = int(np.log2(psi.shape[0]))
    order = copy.deepcopy(ind_A)
    order.extend([p for p in range(num) if p not in ind_A])
    psi = permutebits(psi, order)
    if diag:
        psi = np.abs(psi)**2
        psi = np.reshape(psi, [2**len(ind_A), 2**(num-len(ind_A))])
        psi = np.sum(psi, axis=1)
        return psi        
    else:
        psi = np.reshape(psi, [2**len(ind_A), 2**(num-len(ind_A))]) 
        rho = psi @ np.conj(np.transpose(psi))
        return rho

def py_simulate(qc: QuantumCircuit, 
             state_ini: np.ndarray = np.array([])):
    
    used_qubits = qc.used_qubits
    num = len(used_qubits)
    if not state_ini:
        psi = np.zeros(2**num)
        psi[0] = 1

    else:
        psi = state_ini

    for gate in qc.gates:   
        if not ((isinstance(gate, Delay)) or (isinstance(gate, Barrier))): 
            op = global_op(gate, used_qubits)
            psi = op @ psi


    return psi

def measure_qubits(psi: np.ndarray, qubits: List[int], shots: int): #TODO: Simulate collapsed wave function after measurement
    """Measure qubits with state vector psi
       Args:
            psi: state vector
            qubits: qubits measured
            shots: number of measurments
    """
    num = int(np.log2(psi.shape[0]))
    psi = permutebits(psi, range(num)[::-1])
    probs = ptrace(psi, qubits)

    counts = {}
    for i in range(shots):
        res = np.random.choice(np.arange(len(probs)), p=probs)
        bitstring = bin(res)[2:].zfill(len(qubits))
        if bitstring not in counts.keys():
            counts[bitstring] = 1
        else:
            counts[bitstring] += 1 
    
    return counts

def py_expectation(pypsi : np.ndarray, paulis:list[PauliOp]):
    qubit_num = int(np.log2(len(pypsi)))
    res = []
    for pauli in paulis:
        mat = pauli.get_matrix(qubit_num)
        obs =  np.real(np.conj(pypsi) @ mat @ pypsi.T)
        res.append(obs * pauli.coeff)

    return res