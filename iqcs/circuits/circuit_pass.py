from .quantum_circuit import QuantumCircuit
from .gatelib import *

def decompose_ccx(q1, q2, q3):
    decomposed_gates = []
    decomposed_gates.append(H(q3))
    decomposed_gates.append(CX(q2, q3))
    decomposed_gates.append(Tdg(q3))
    decomposed_gates.append(CX(q1, q3))
    decomposed_gates.append(T(q3))
    decomposed_gates.append(CX(q2, q3))
    decomposed_gates.append(Tdg(q3))
    decomposed_gates.append(CX(q1, q3))
    decomposed_gates.append(T(q2))
    decomposed_gates.append(T(q3))
    decomposed_gates.append(H(q3))
    decomposed_gates.append(CX(q1, q2))
    decomposed_gates.append(T(q1))
    decomposed_gates.append(Tdg(q2))
    decomposed_gates.append(CX(q1, q2))

    return decomposed_gates

def decompose_cswap(q1, q2, q3):
    decomposed_gates = []
    decomposed_gates.append(CX(q3, q2))
    decomposed_gates += decompose_ccx(q1, q2, q3)
    decomposed_gates.append(CX(q3, q2))
    return decomposed_gates

def decompose_3bit_gate(circuit: QuantumCircuit):
    decomposed_gates = []
    for gate in circuit.operations:
        if gate.name == "CSWAP":
            q1, q2, q3 = gate.pos
            decomposed_gates += decompose_cswap(q1, q2, q3)
        elif gate.name == "CCX":
            q1, q2, q3 = gate.pos
            decomposed_gates += decompose_ccx(q1, q2, q3)
        else:
            decomposed_gates.append(gate)

    circuit.operations = decomposed_gates
