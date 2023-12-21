
from iqcs import QuantumCircuit
from typing import Iterable
from math import pi
from iqcs.circuits.gatelib import GATE_LIB, SX, SY
from copy import deepcopy

def json_to_circuit(gateList, n_qubits):
    qc = QuantumCircuit(n_qubits)
    measures = []
    cbits = []
    for item in gateList:
        args = deepcopy(item['qubits'])
        if 'angle' in item.keys():
            args += item['angle']

        if  item['name']=='measure':
            measures.append(item['qubits'][0])
            cbits.append(item['memory'][0])
        elif item['name']=='x_half':
            qc << (SX(item['qubits'][0]))
        elif item['name']=='y_half':
             qc << (SY(item['qubits'][0]))
        else:
            qc << (GATE_LIB[item['name']](*args))

    qc.measure(measures, cbits)
    return qc

def circuit_to_json(qc):
    gate_list = []
    for gate in qc.operations:
        gate_name = gate.name.lower()
        if gate.name == "sx":
            gate_name = "x_half"
        if gate.name == "sy":
            gate_name = "y_half"

        if gate_name not in ["delay", "barrier"]:
            gate_dict = {"name"  : gate_name}
            gate_dict["qubits"] = deepcopy(gate.pos)
          
            if not gate.paras == None:
                gate_dict["angle"] = deepcopy(gate.paras)
               
            
            gate_list.append(gate_dict)
    
    for qbit, cbit in qc.measures.items():
        gate_list.append({"memory":[cbit], "name":"measure", "qubits":[qbit]})
        
    return gate_list