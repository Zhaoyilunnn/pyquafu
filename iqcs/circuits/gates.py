from typing import List ,Union, Callable
import numpy as np
from .utils import reorder_matrix, extract_float
from .operations import Operation
import copy

from .parameters import Parameter, ParameterExpression
ParameterType=Union[float, Parameter, ParameterExpression]


class QuantumGate(Operation):
    def __init__(self, 
                 name: str, 
                 pos: List[int], 
                 paras: List[ParameterType], 
                 matrix:Union[np.ndarray,Callable]):
        
        self._name = name
        self._pos = pos
        self.paras = paras

        self._raw_matrix = matrix

    @property
    def name(self):
        return self._name 
    
    @name.setter
    def name(self, __name):
        self._name = __name
    
    @property
    def pos(self):
        return  self._pos
    
    @pos.setter
    def pos(self, __pos):
        self._pos = copy.deepcopy(__pos)
    
    @property
    def _paras(self):
        return extract_float(self.paras)
    
    
    @property
    def matrix(self):
        raw_mat = self._raw_matrix
        if isinstance(self._raw_matrix, Callable):
            raw_mat = self._raw_matrix(self._paras)

        if len(self.pos) > 1:
            return reorder_matrix(raw_mat, self.pos)
        else:
            return raw_mat

    @property
    def symbol(self):
        #TODO: Use latex repr for Parameter
        if len(self.paras) > 0:
            symbol = "%s(" %self.name + ",".join(["%.3f" %para for para in self._paras]) + ")"
            return symbol
        else:
            return "%s" %self.name

    @symbol.setter
    def symbol(self, symbol):
        self.symbol = symbol

    def __str__(self):
        properties_names = ['pos', 'paras']
        properties_values = [getattr(self, x) for x in properties_names]
        return "%s:\n%s" % (self.__class__.__name__, '\n'.join(
            [f"{x} = {repr(properties_values[i])}" for i, x in enumerate(properties_names)]))

    def __repr__(self):
        return f"{self.__class__.__name__}"

    def to_qasm(self):
        #TODO:qasm parser support variables
        qstr = "%s" %self.name.lower()
        if self.paras:
            qstr += "(" + ",".join(["%s" %para for para in self._paras]) + ")"
        qstr += " "
        qstr += ",".join(["q[%d]" % p for p in self.pos])
        return qstr
    
    def power(self, n) -> "QuantumGate":
        name = self.name
        if self.name in ["RX", "RY", "RZ", "P", "RZZ", "RXX", "RYY"]:
            return QuantumGate(name, self.pos, [n * self.paras[0]], self._raw_matrix)
        else:
            name = self.name + "^%d" %n
            raw_matrix = self._raw_matrix 
            if isinstance(self._raw_matrix, Callable):
                raw_matrix = lambda paras: np.linalg.matrix_power(self._raw_matrix(paras), n)
            else:
                raw_matrix = np.linalg.matrix_power(self._raw_matrix, n)
            return QuantumGate(name, self.pos, self.paras, raw_matrix)
        
    def dagger(self) -> "QuantumGate":
        name = self.name
        if name in ["RX", "RY", "RZ", "P", "RZZ", "RXX", "RYY"]:
            return QuantumGate(name, self.pos, [-self.paras[0]], self._raw_matrix)
        else:
            name = self.name + "^†"
            raw_matrix = self._raw_matrix 
            if isinstance(self._raw_matrix, Callable):
                raw_matrix = lambda paras: self._raw_matrix(paras).conj().T
            else:
                raw_matrix = raw_matrix.conj().T
            return QuantumGate(name, self.pos, self.paras, raw_matrix)

    def add_controls(self, ctrls) -> "QuantumGate":
        from .gatelib import GATE_LIB
        if self.name in ["X" ,"Y" ,"Z", "RX", "RY", "RZ"]:
            args = self.pos + self.paras
            cname = "mc" + self.name.lower()
            cop = GATE_LIB[cname](ctrls, *args)
            return cop
        else:
            return ControlledU("C"+self.name, ctrls, self)
        
    def _get_raw_matrix(self, reverse_order=False):
        raw_mat = self._raw_matrix
        if isinstance(self._raw_matrix, Callable):
            raw_mat = self._raw_matrix(self._paras)
        if reverse_order and len(self.pos) > 1:
            return reorder_matrix(raw_mat, np.arange(len(self.pos))[::-1])
        else:
           return raw_mat
    

class ControlledGate(QuantumGate):
    """ Controlled gate class, where the matrix act non-trivaly on target qubits"""
    def __init__(self, name:str, 
                 targ_name:str, 
                 ctrls: List[int], 
                 targs: List[int], 
                 paras : List[float], 
                 targ_matrix:Union[np.ndarray,Callable]):
        self.ctrls = copy.deepcopy(ctrls)
        self.targs = copy.deepcopy(targs)
        self._targ_name = targ_name
        super().__init__(name, ctrls+targs, paras, targ_matrix)
        self._targ_matrix =  targ_matrix
        self._raw_matrix = self._rawmatfunc 
        
    @property
    def symbol(self):
        if len(self.paras) > 0:
            symbol = "%s(" %self._targ_name + ",".join(["%.3f" %para for para in self._paras]) + ")"
            return symbol
        else:
            return "%s" %self._targ_name
    
    def _rawmatfunc(self, paras:List[float]):
        targ_dim = 2**(len(self.targs))
        qnum = len(self.pos)
        dim = 2**(qnum)
        raw_matrix =  np.zeros((dim , dim), dtype=complex)
        targ_matrix = self._targ_matrix
        if isinstance(self._targ_matrix, Callable):
            targ_matrix = self._targ_matrix(paras)

        if targ_matrix.shape[0] != targ_dim:
            raise ValueError("Dimension dismatch")
        else:
            control_dim = 2**len(self.pos) - targ_dim
            for i in range(control_dim):
                raw_matrix[i, i] = 1.
            
            raw_matrix[control_dim:, control_dim:] = targ_matrix

        return raw_matrix
    
    def power(self, n) -> "ControlledGate":
        name = self.name
        if self._targ_name in ["RX", "RY", "RZ", "P", "RZZ", "RXX", "RYY"]:
            return ControlledGate(name, self._targ_name, self.ctrls, self.targs, [self.paras[0] * n], self._targ_matrix)
        else:
            name = self.name + "^%d" %n
            targ_matrix = self._targ_matrix 
            if isinstance(self._targ_matrix, Callable):
                targ_matrix = lambda paras: np.linalg.matrix_power(self._targ_matrix(paras), n)
            else:
                targ_matrix = np.linalg.matrix_power(self._targ_matrix, n)
            return ControlledGate(name, self._targ_name+"^%d" %n, self.ctrls, self.targs, self.paras, targ_matrix)
    
        
    def dagger(self) -> "ControlledGate":
        name = self.name
        if self._targ_name in ["RX", "RY", "RZ", "P", "RZZ", "RXX", "RYY"]:
            return ControlledGate(name, self._targ_name, self.ctrls, self.targs, [-self.paras[0]], self._targ_matrix)
        else:
            name = self.name + "^†"
            targ_matrix = self._targ_matrix 
            if isinstance(self._targ_matrix, Callable):
                targ_matrix = lambda paras: self._targ_matrix(paras).conj().T
            else:
                targ_matrix = targ_matrix.conj().T

            return ControlledGate(name, self._targ_name+"^+", self.ctrls, self.targs, self.paras, targ_matrix)
    
    def add_controls(self, ctrls) -> "ControlledGate":
        from .gatelib import GATE_LIB
        if self.name in ["CX", "CY", "CZ", "CRX", "CRY", "CRZ"]:
            args = self.targs + self.paras
            cname = "mc" + self._targ_name.lower()
            cop = GATE_LIB[cname](ctrls+self.ctrls, *args)
            return cop
        else:
            cop = ControlledGate("MC"+self._targ_name, self._targ_name, ctrls+self.ctrls, self.targs, self.paras, self._targ_matrix)
            return cop
        
    def _get_targ_matrix(self, reverse_order=False):
        targ_mat = self._targ_matrix
        if isinstance(self._targ_matrix, Callable):
            targ_mat = self._targ_matrix(self._paras)
        if reverse_order and (len(self.targs) > 1): 
            return reorder_matrix(targ_mat, np.array(range(len(self.targs))[::-1]))
        else:
            return targ_mat
    
class ControlledU(ControlledGate):
    def __init__(self, name, ctrls: List[int], U: QuantumGate):
        self.targ_gate = U
        targs = U.pos
        super().__init__(name, U.name, ctrls, targs, U.paras, targ_matrix=self.targ_gate._raw_matrix)

class CircuitWrapper(QuantumGate):
    def __init__(self, name:str, circ, qbits=[]):
        self.name  = name
        self.pos = list(range(circ.num))
        self.circuit = copy.deepcopy(circ)

#TODO:Handle wrapper paras
        # if hasattr(circ, "paras"):
        #     self._paras = circ.paras
        # else:
        #     self._paras = []
        #     for op in self.circuit.operations:
        #         self._paras.extend(op.paras)
        
        if qbits:
            self._reallocate(qbits)

    # @property
    # def paras(self):
    #     return self._paras
    
    # @paras.setter
    # def paras(self, __paras):
    #     self._paras = __paras
    #     self.circuit.paras = __paras
    
    def _reallocate(self, qbits):
        num = max(self.circuit.num-1, max(qbits))+1
        self.pos = qbits
        self.circuit._reallocate(num, qbits)
    
    @property
    def symbol(self):
        return "%s" %self.name

    def add_controls(self, ctrls:List[int]=[]) -> QuantumGate:
        return ControlCircuitWrapper("MC"+self.name, self, ctrls)

    def power(self, n:int):
        self.name +="^%d" %n
        self.circuit = self.circuit.power(n)
        return self

    def dagger(self):
        self.name += "^†"
        self.circuit = self.circuit.dagger()
        return self

    def to_qasm(self):
        qasm = ''
        for operation in self.circuit.operations:
            qasm += operation.to_qasm() + ";\n"
        return qasm
    
class ControlCircuitWrapper(CircuitWrapper):
    def __init__(self, name:str, circwrp:CircuitWrapper, ctrls:List[int]):
        self.name = name
        self.ctrls = ctrls
        self.targs = circwrp.pos
        self.circuit = circwrp.circuit.add_controls(len(ctrls), ctrls, self.targs)
        self.pos = list(range(self.circuit.num))
        self._targ_name = circwrp.name

    
    @property
    def symbol(self):
        return "%s" %self._targ_name

    def power(self, n: int):
        self._targ_name += "^%d" %n
        return super().power(n)

    def dagger(self):
        self.name += "^†"
        return super().dagger()
    
    def _reallocate(self, qbits):
        num = max(self.circuit.num-1, max(qbits))+1
        self.pos = qbits
        self.circuit._reallocate(num, qbits)
        qbits_map = dict(zip(range(len(qbits)), qbits))
        for i in range(len(self.ctrls)):
            self.ctrls[i] = qbits_map[self.ctrls[i]]
        
        for i in range(len(self.targs)):
            self.targs[i] = qbits_map[self.targs[i]]
