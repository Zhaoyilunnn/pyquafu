
from typing import Union

from .iqsimu import simulate_statevec, simulate_mps, expect_statevec, applyop_statevec
from iqcs import QuantumCircuit
from .simuresults import SimuResult
import numpy as np
from ..exceptions import IqcsError
from ..circuits.circuit_pass import decompose_3bit_gate
from ..utils.paulis import Hamiltonian
from typing import List
from abc  import ABC, abstractmethod
from ..circuits.gatelib import GATE_LIB
from ..circuits.gates import QuantumGate, CircuitWrapper

class Simulator(ABC):
    @abstractmethod
    def run(self):
        raise NotImplementedError

class SVSimulator(Simulator):
    def __init__(self, use_gpu:bool=False, use_custatevec:bool=False, merge_gates=True, max_merge_size=5):
        self.use_gpu  = use_gpu
        self.use_custatevec = use_custatevec
        self.merge_gates = merge_gates
        self.max_merge_size = max_merge_size
    
    def config(self, **kwargs):
        for attr, value in kwargs.items():
            if hasattr(self, attr):
                setattr(self, attr, value)
            else:
                raise ValueError("No such attribute")

    def _apply_op(self, op, psi):
        #TODO: support GPU
        if isinstance(op, CircuitWrapper):
            psi = self.run(op.circuit, psi)["statevector"]
            return psi
        elif isinstance(op, QuantumGate):
            psi = applyop_statevec(op, psi)
            return psi
        else:
            raise NotImplementedError
    
    def _apply_hamil(self, hamil, psi):
        psi_out = np.zeros(len(psi), dtype=complex) 
        for pauli in hamil.paulis:
            psi1 = np.copy(psi)
            for name, pos in zip(pauli.paulistr, pauli.pos):
                op = GATE_LIB[name.lower()](pos)
                psi1 = self._apply_op(op, psi1)
            psi1 = psi1 * pauli.coeff
            psi_out += psi1
    
        return psi_out

    def run(self, qc : QuantumCircuit, psi : np.ndarray= np.array([]), shots:int=0, hamiltonian:Hamiltonian=None):
        res_dict = {}
        num  = qc.num
        if num <= 20:
            self.merge_gates = False

        if self.use_gpu:
            if self.use_custatevec:
                try:
                    from .iqsimu import simulate_statevec_custate
                except ImportError:
                    raise IqcsError(" pyiqcs is installed without cuquantum support")
                psi = simulate_statevec_custate(qc, psi)
            else:
                try:
                    from .iqsimu import simulate_statevec_gpu
                except ImportError:
                    raise IqcsError("You are not using the GPU version of pyiqcs")
                psi = simulate_statevec_gpu(qc, psi)
        else:
            res_dict = simulate_statevec(qc, shots, psi, self.merge_gates, self.max_merge_size)
            
            if hamiltonian:
                paulis = hamiltonian.paulis
                res = expect_statevec(res_dict["statevector"], paulis)
                for i in range(len(paulis)):
                    res[i] *= paulis[i].coeff
                res_dict["pauli_expects"] = res
            else:
                res_dict["pauli_expects"] = []
            res_dict["measures"] =  qc.measures
            return SimuResult(res_dict)
        
class MPSSimulator(Simulator):
    def  __init__(self,  cutoff=1e-16, save_statevec=False):
        self.cutoff = cutoff
        self.save_statevec = save_statevec
    
    def config(self, **kwargs):
        for attr, value in kwargs.items():
            if hasattr(self, attr):
                setattr(self, attr, value)
            else:
                raise ValueError("No such attribute")

    def run(self, qc : QuantumCircuit, shots:int=0, hamiltonian:Hamiltonian=None):
        decompose_3bit_gate(qc)
        paulis = []
        if hamiltonian:
            paulis = hamiltonian.paulis

        res_dict =  simulate_mps(qc, cutoff=self.cutoff, shots=shots, save_statevec=self.save_statevec, paulis=paulis)

        for i in range(len(paulis)):
            res_dict["pauli_expects"][i] *= paulis[i].coeff
        res_dict["measures"] =  qc.measures 
        return SimuResult(res_dict)


