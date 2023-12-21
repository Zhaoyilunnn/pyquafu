
from ..gates import QuantumGate, ControlledGate
from .gate_mats import *

from ..gates import ParameterType

GATE_LIB = {}

def register_gate(cls):
    cls_name = cls.__name__.lower()
    def register(cls):
        GATE_LIB[cls_name] = cls
        return cls
    return register(cls)

@register_gate
class ID(QuantumGate):
    name = "ID"
    _raw_matrix = ID_MAT
    paras = [] 
    def __init__(self, pos: int):
       self.pos =  [pos]

@register_gate
class H(QuantumGate):
    name = "H"
    _raw_matrix = H_MAT
    paras = [] 
    def __init__(self, pos: int):
       self.pos =  [pos]

@register_gate
class X(QuantumGate):
    name = "X"
    _raw_matrix = X_MAT
    paras = [] 
    def __init__(self, pos: int):
       self.pos =  [pos]

@register_gate
class Y(QuantumGate):
    name = "Y"
    _raw_matrix = Y_MAT
    paras = [] 
    def __init__(self, pos: int):
       self.pos =  [pos]

@register_gate
class Z(QuantumGate):
    name = "Z"
    _raw_matrix = Z_MAT
    paras = [] 
    def __init__(self, pos: int):
       self.pos =  [pos]

@register_gate
class S(QuantumGate):
    name = "S"
    _raw_matrix = S_MAT
    paras = [] 
    def __init__(self, pos: int):
       self.pos =  [pos]

@register_gate
class Sdg(QuantumGate):
    name = "Sdg"
    _raw_matrix = S_MAT.conj().T
    paras = [] 
    def __init__(self, pos: int):
       self.pos =  [pos]

@register_gate
class T(QuantumGate):
    name = "T"
    _raw_matrix = T_MAT
    paras = [] 
    def __init__(self, pos: int):
       self.pos =  [pos]

@register_gate
class Tdg(QuantumGate):
    name = "Tdg"
    _raw_matrix = T_MAT.conj().T
    paras = [] 
    def __init__(self, pos: int):
       self.pos =  [pos]

@register_gate
class W(QuantumGate):
    name = "W"
    _raw_matrix = W_MAT
    paras = [] 
    def __init__(self, pos: int):
       self.pos =  [pos]
    
    def to_qasm(self):
        q = self.pos[0]
        return "rz(-pi/4) q[%d];\nrx(pi) q[%d];\nrz(pi/4) q[%d]"  %(q, q, q)

@register_gate
class SX(QuantumGate):
    name = "SX"
    _raw_matrix = SX_MAT
    paras = [] 
    def __init__(self, pos: int):
       self.pos =  [pos]

@register_gate
class SY(QuantumGate):
    name = "SY"
    _raw_matrix = SY_MAT
    paras = [] 
    def __init__(self, pos: int):
       self.pos =  [pos]

    def to_qasm(self):
        return "ry(pi/2) q[%d]" %(self.pos[0])

@register_gate
class SW(QuantumGate):
    name = "SW"
    _raw_matrix = SW_MAT
    paras = [] 
    def __init__(self, pos: int):
       self.pos =  [pos]

    def to_qasm(self):
        q = self.pos[0]
        return "rz(-pi/4) q[%d];\nrx(pi/2) q[%d];\nrz(pi/4) q[%d]"  %(q, q, q)


@register_gate
class RX(QuantumGate):
    def __init__(self, pos: int, theta:ParameterType):
        super().__init__("RX", [pos], [theta], rxmatrix)

@register_gate
class RY(QuantumGate):
    def __init__(self, pos: int, theta:ParameterType):
        super().__init__("RY", [pos], [theta], rymatrix)

@register_gate
class RZ(QuantumGate):
    def __init__(self, pos: int, theta:ParameterType):
        super().__init__("RZ", [pos], [theta], rzmatrix)

# @register_gate
# class U2(QuantumGate):
#     def __init__(self, pos: int, phi: float, _lambda: float):
#         super().__init__("U2", [pos], [phi, _lambda], u2matrix(phi, _lambda))

@register_gate
class U3(QuantumGate):
    def __init__(self, pos: int, theta : ParameterType, phi: ParameterType, _lambda: ParameterType):
        super().__init__("U3", [pos], [theta, phi, _lambda], matrix = u3matrix)

@register_gate
class P(QuantumGate):
    def __init__(self, pos: int, _lambda:ParameterType):
        super().__init__("P", [pos], [_lambda], pmatrix)

@register_gate
class ISWAP(QuantumGate):
    name = "ISWAP"
    matrix = ISWAP_MAT
    _raw_matrix = ISWAP_MAT
    paras = [] 
    def __init__(self, q1:int, q2:int):
       self.pos =  [q1,  q2]

@register_gate
class SWAP(QuantumGate):
    name = "SWAP"
    matrix = SWAP_MAT
    _raw_matrix = SWAP_MAT
    paras = [] 
    def __init__(self, q1:int, q2:int):
       self.pos =  [q1,  q2]
    
    @property
    def symbol(self):
        return "x"

@register_gate
class RXX(QuantumGate):
    def __init__(self, q1:int, q2:int, theta:ParameterType):
        super().__init__("RXX", [q1, q2], [theta], rxxmatrix)
    
    
@register_gate
class RYY(QuantumGate):
    def __init__(self, q1:int, q2:int, theta:ParameterType):
        super().__init__("RYY", [q1, q2], [theta], ryymatrix)

    
@register_gate
class RZZ(QuantumGate):
    def __init__(self, q1:int, q2:int, theta:ParameterType):
        super().__init__("RZZ", [q1, q2], [theta], rzzmatrix)

  
@register_gate
class CX(ControlledGate):
    name  = "CX"
    _targ_name = "X"
    _targ_matrix = X_MAT
    _raw_matrix = CX_MAT
    paras = []
    def __init__(self, ctrl:int, targ:int):
        assert ctrl != targ
        self.ctrls  = [ctrl]
        self.targs = [targ]
        self.pos = self.ctrls + self.targs

    @property
    def symbol(self):
        return "+"

@register_gate
class CY(ControlledGate):
    name  = "CY"
    _targ_name = "Y"
    _targ_matrix = Y_MAT
    _raw_matrix = CY_MAT
    paras = []
    def __init__(self, ctrl:int, targ:int):
        assert ctrl != targ
        self.ctrls  = [ctrl]
        self.targs = [targ]
        self.pos = self.ctrls + self.targs

@register_gate
class CZ(ControlledGate):
    name  = "CZ"
    _targ_name = "Z"
    _targ_matrix = Z_MAT
    _raw_matrix = CZ_MAT
    paras = []
    def __init__(self, ctrl:int, targ:int):
        assert ctrl != targ
        self.ctrls  = [ctrl]
        self.targs = [targ]
        self.pos = self.ctrls + self.targs

@register_gate
class CS(ControlledGate):
    name  = "CS"
    _targ_name = "S"
    _targ_matrix = S_MAT
    _raw_matrix = CS_MAT
    paras = []
    def __init__(self, ctrl:int, targ:int):
        assert ctrl != targ
        self.ctrls  = [ctrl]
        self.targs = [targ]
        self.pos = self.ctrls + self.targs

    def to_qasm(self):
        return "cp(pi/2) " + "q[%d],q[%d]" % (self.pos[0], self.pos[1])
    
@register_gate
class CT(ControlledGate):
    name  = "CT"
    _targ_name = "T"
    _targ_matrix = T_MAT
    _raw_matrix = CT_MAT
    paras = []
    def __init__(self, ctrl:int, targ:int):
        assert ctrl != targ
        self.ctrls  = [ctrl]
        self.targs = [targ]
        self.pos = self.ctrls + self.targs

    def to_qasm(self):
        return "cp(pi/4) " + "q[%d],q[%d]" % (self.pos[0], self.pos[1])

@register_gate
class CP(ControlledGate):
    def __init__(self, ctrl:int, targ:int, _lambda:ParameterType):
        super().__init__("CP", "P", [ctrl], [targ], [_lambda], pmatrix)

@register_gate
class CRX(ControlledGate):
    def __init__(self, ctrl:int, targ:int, theta:ParameterType):
        super().__init__("CRX", "RX", [ctrl], [targ], [theta], rxmatrix)

@register_gate
class CRY(ControlledGate):
    def __init__(self, ctrl:int, targ:int, theta:ParameterType):
        super().__init__("CRY", "RY", [ctrl], [targ], [theta], rymatrix)

@register_gate
class CRZ(ControlledGate):
    def __init__(self, ctrl:int, targ:int, theta:ParameterType):
        super().__init__("CRZ", "RZ", [ctrl], [targ], [theta], rzmatrix)

@register_gate
class CCX(ControlledGate):
    def __init__(self, ctrl1:int, ctrl2:int, targ:int):
        super().__init__("CCX", "X", [ctrl1, ctrl2], [targ], [], X_MAT)

@register_gate
class CSWAP(ControlledGate):
    def __init__(self, ctrl:int, targ1:int, targ2:int):
        super().__init__("CSWAP", "SWAP", [ctrl], [targ1, targ2], [], SWAP_MAT)

@register_gate
class MCX(ControlledGate):
    def __init__(self, ctrls:List[int], targ:int):
        super().__init__("MCX", "X", ctrls, [targ], [], X_MAT)
    
    @property
    def symbol(self):
        return "+"

@register_gate
class MCY(ControlledGate):
    def __init__(self, ctrls:List[int], targ:int):
        super().__init__("MCY", "Y", ctrls, [targ], [], Y_MAT)

@register_gate
class MCZ(ControlledGate):
    def __init__(self, ctrls:List[int], targ:int):
        super().__init__("MCZ", "Z", ctrls, [targ], [], Z_MAT)

@register_gate
class MCRX(ControlledGate):
    def __init__(self, ctrls:List[int], targ:int, theta:ParameterType):
        super().__init__("MCRX", "RX", ctrls, [targ], [theta], rxmatrix)

@register_gate
class MCRY(ControlledGate):
    def __init__(self, ctrls:List[int], targ:int, theta:ParameterType):
        super().__init__("MCRY", "RY", ctrls, [targ], [theta], rymatrix)

@register_gate
class MCRZ(ControlledGate):
    def __init__(self, ctrls:List[int], targ:int, theta:ParameterType):
        super().__init__("MCRZ", "RZ", ctrls, [targ], [theta], rzmatrix)
 

GATE_LIB["cnot"] = CX
GATE_LIB["toffoli"] = CCX
GATE_LIB["fredkin"] = CSWAP