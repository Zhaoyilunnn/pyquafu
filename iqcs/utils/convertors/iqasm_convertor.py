
from iqcs.iqasm.iqasm_parser import IqasmParser, QregNode
from iqcs.iqasm.instruction_node import InstructionNode, node_to_gate
from iqcs import QuantumCircuit

def qasmfile_to_circuit(qasmfile):
    qasmstr = ""
    with open(qasmfile) as f:
        for line in f.readlines():
            qasmstr += line

    return qasm_to_circuit(qasmstr)
 
def qasm_to_circuit(qasm):
    """TODO:Support u1, u2, u3 gate"""
    parser = IqasmParser()     
    nodes = parser.parse(qasm)

    n = 0
    operations = []
    measures = {}
    for node in nodes:
        if isinstance(node, QregNode):
            n = node.n
        if isinstance(node, InstructionNode):
            if node.name == "measure":
                for q, c in zip(node.pos.keys(), node.pos.values()):
                    measures[q] = c
            else:
                operations.append(node_to_gate(node))

    q = QuantumCircuit(n)
    q.operations = operations
    q.measures = measures
    return q

