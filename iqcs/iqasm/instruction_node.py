from typing import Dict, Any, List, Union
import dataclasses 
from iqcs.circuits.gatelib import GATE_LIB

@dataclasses.dataclass
class InstructionNode: 
    name:Any            # gate.name
    pos:Union[List[Any], Dict[Any,Any]]       # gate.pos |  Dict[Any,Any] for measure
    paras:List[Any]        # gate.paras
    # matrix:List[Any]   # for gate in [QuantumGate]
    duration:int       # for gate in Delay in iqcs
    unit:str           # for gate in Delay in iqcs
    label:str
        
    def __hash__(self):
        return hash((type(self.name), tuple(self.pos) ,self.label))
    
    def __str__(self):
        if self.name == 'measure':
            args = ','.join(str(q) for q in self.pos.keys())
            args += f'=>{",".join(str(c) for c in self.pos.values())}'
        else:    
            args = ','.join(str(q) for q in self.pos)
            
        if self.paras == None:
            return f'{self.label}{{{self.name}({args})}}' # no paras
        else:
            # if self.paras not a list, then make it a list  of str of .3f float
            if not isinstance(self.paras, list):
                formatted_paras = [f'{self.paras:.3f}']
            else:
                formatted_paras = [f'{p:.3f}' for p in self.paras]  
                
            formatted_paras_str = ','.join(formatted_paras)
            
            return f'{self.label}{{{self.name}({args})}}({formatted_paras_str})'
    
    def __repr__(self): 
        return str(self)
    
def node_to_gate(gate_in_dag):
    """
    transform gate in dag graph, to gate in circuit which can be added to circuit

    Args:
        gate_in_dag: a node in dag graph , gate_in_dag is a GateWrapper object. 
            in GateWrapper, gate_in_dag.name is uppercase, gate_in_dag.pos is a list or a dict
            gate_transform support gate with one qubit or more qubits, not measures!
            and you should exculde nodes [-1 ,float('inf') , measure_gate] in dag graph

    Returns:
        gate: gate  which can be added to circuit in iqcs
    """

    gate_name = gate_in_dag.name.lower()
    gate_class = GATE_LIB.get(gate_name)
    if not gate_class:
        raise ValueError("gate %s is not supported" %gate_name)

    if gate_name == "barrier":
        return gate_class(gate_in_dag.pos)

    # 从gate_in_dag获取参数列表
    args = gate_in_dag.pos
    if gate_in_dag.paras:
        args += gate_in_dag.paras

    # 处理 gate.duration 和 gate.unit
    if gate_name in ["delay"]:
        args.append(gate_in_dag.duration)
        args.append(gate_in_dag.unit)

    # 处理多量子比特门
    if gate_name in ["mcx", "mcy", "mcz"]:
        control_qubits = gate_in_dag.pos[:-1]
        target_qubit = gate_in_dag.pos[-1]
        return gate_class(control_qubits, target_qubit)

    # measure node
    if gate_name == "measure":
        return gate_class(gate_in_dag.pos)
    
    return gate_class(*args)