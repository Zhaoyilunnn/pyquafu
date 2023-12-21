from typing import Iterable, List, Union
import numpy as np
from .gates import QuantumGate, ControlledGate, ControlledU, CircuitWrapper, ControlCircuitWrapper
from .gatelib import GATE_LIB
from .operations import Operation, Barrier, Delay 
from ..exceptions import CircuitError
import copy
from .parameters import Parameter, ParameterExpression

class QuantumCircuit(object):
    def __init__(self, num: int, name=""):
        """
        Initialize a QuantumCircuit object
        
        Args:   
            num (int): Total qubit number used
        """
        self.num = num
        self.operations = []
        self.openqasm = ""
        self.circuit = []
        self.measures = dict(zip(range(num), range(num)))
        self._used_qubits = []
        self.name = name
        self._parameter_grads = {}
        self._variables = []
        self._has_wrap = False
    @property
    def used_qubits(self) -> List:
        self.layered_circuit()
        return self._used_qubits

    def __lshift__(self, operation :Operation):
        if max(operation.pos) >= self.num:
            raise CircuitError("Operation act on qubit that not allocated")
        self.operations.append(operation)
        if isinstance(operation, CircuitWrapper) or isinstance(operation, ControlCircuitWrapper):
            self._has_wrap = True
        return self
    
    def get_parameter_grads(self):
        if self._has_wrap:
            print("warning: The circuit has wrapped gates, it will unwarp automaticllay")
            self.unwrap()
        self._parameter_grads = {}
        for i, op in enumerate(self.operations):
            for j, para in enumerate(op.paras):
                if isinstance(para, Parameter):
                    if para not in self._parameter_grads.keys():
                       self._parameter_grads[para] = [[(i, j), 1.]]
                    else:
                       self._parameter_grads[para].append([(i, j), 1.])

                elif isinstance(para, ParameterExpression):
                    para_grads = para.grad()
                    for var in para._variables:
                        if var not in self._parameter_grads.keys():
                            self._parameter_grads[var] = [[(i,j),  para_grads[para._variables[var]]]]
                        else:
                            self._parameter_grads[var].append([(i,j), para_grads[para._variables[var]]])
        self._variables = list(self._parameter_grads.keys())
        return self._parameter_grads
    
    def _calc_parameter_grads(self):
        """
        Update parameter gradient value, not variables
        """
        for var, grads in self._parameter_grads.items():
             for item in grads:
                i, j = item[0]
                para = self.operations[i].paras[j]
                if isinstance(para, ParameterExpression):
                    para_grads = para.grad()
                    item[1] = para_grads[para._variables[var]]

        return self._parameter_grads
    
    @property
    def variables(self):
        self._variables = list(self.get_parameter_grads().keys())
        return self._variables
    
    def update_params(self, values):
        """
        Update variables' value, not variables
        """

        for i in range(len(values)):
            self._variables[i].value = values[i]
    
    def _reallocate(self, num, qbits:List[int]):
        assert self.num == len(qbits)
        if max(qbits) > num:
            raise CircuitError("Bad allocation")
        
        self.num = num
        qbits_map = dict(zip(range(len(qbits)), qbits))
        operations = self.operations
        for op in operations:
            for i in range(len(op.pos)):
                op.pos[i] = qbits_map[op.pos[i]]
            
            if isinstance(op, ControlledGate):
                for i in range(len(op.ctrls)):
                    op.ctrls[i] = qbits_map[op.ctrls[i]]
                
                for i in range(len(op.targs)):
                    op.targs[i] = qbits_map[op.targs[i]]
                
                if isinstance(op, ControlledU):
                    for i in range(len(op.targ_gate.pos)):
                        op.targ_gate.pos[i] = qbits_map[op.targ_gate.pos[i]]


    def add_controls(self, ctrlnum, ctrls:List[int]=[], targs:List[int]=[], inplace=False)->"QuantumCircuit":
        num = 0
        operations = []
        if len(ctrls+targs) == 0:
            ctrls = list(np.arange(ctrlnum) + self.num)
            num = self.num + ctrlnum
            operations = self.operations
        else:
            if len(ctrls) == 0 or len(targs) == 0:
                raise ValueError("Must provide both ctrls and targs")
            else:
                assert len(targs) == self.num
                assert len(ctrls) == ctrlnum  
                num = max(ctrls+targs)+1 
                if inplace:       
                    self._reallocate(num, targs)
                else:
                    temp = copy.deepcopy(self)
                    temp._reallocate(num, targs)
                    operations = temp.operations
        
        if inplace:
            for op in self.operations:
               op = op.add_controls(ctrls)
            return self
        else:
            qc = QuantumCircuit(num)
            for op in operations:
                qc << op.add_controls(ctrls)
            return qc

    def join(self, 
             qc : ["QuantumCircuit", CircuitWrapper, ControlCircuitWrapper], 
             qbits:List[int]=[], 
             inplace=True)->"QuantumCircuit":
        
        num = self.num
        rnum = 0
        if isinstance(qc, QuantumCircuit):
            rnum = qc.num
        else:
            rnum = qc.circuit.num

        if len(qbits) == 0:
            num = max(self.num, rnum)
        else:
            num = max(self.num-1, max(qbits))+1

        nqc = QuantumCircuit(num)
        if inplace:
            self.num = num
            nqc = self
        else:
            for op in self.operations:
                nqc << op
        
        if isinstance(qc, QuantumCircuit):
            if qbits:
                newqc = copy.deepcopy(qc)
                newqc._reallocate(num, qbits)
                for op in newqc.operations:
                    nqc << op
            else:
                for op in qc.operations:
                    nqc << op
        else:
            if qbits:
                qc._reallocate(qbits)
            nqc << qc
   
        return nqc
       

    def power(self, n, inplace=False):
        if inplace:
            for _ in range(n-1):
                self.operations += self.operation
            return self
        else:
            nq = QuantumCircuit(self.num)
            for _ in range(n):
                for op in self.operations:
                    nq << op
            nq.measures = self.measures
            return nq
        
    def dagger(self, inplace=False):
        if inplace:
            self.operations = self.operations[::-1]
            for op in self.operations:
                op = op.dagger()
            return self
        else:
            nq = QuantumCircuit(self.num)
            for op in self.operations[::-1]:
                nq << op.dagger()
            nq.measures = self.measures
            return nq
        
    def wrap(self, qbits=[]):
        name = self.name if self.name else "Oracle"
        return CircuitWrapper(name, self, qbits)
    
    def unwrap(self):
        operations = []
        for op in self.operations:
            if isinstance(op, CircuitWrapper) or isinstance(op, ControlCircuitWrapper):
                circ = op.circuit.unwrap()
                for op_ in circ.operations:
                    operations.append(op_)
            else:
                operations.append(op)
        self.operations = operations
        self._has_wrap = False
        return self
    
    def layered_circuit(self) -> np.ndarray:
        """
        Make layered circuit from the operation sequence self.operations.

        Returns: 
            A layered list with left justed circuit.
        """
        num = self.num
        operationlist = self.operations
        operationQlist = [[] for i in range(num)]
        used_qubits = []
        for operation in operationlist:
            if isinstance(operation, Delay):
                operationQlist[operation.pos].append(operation)
                if operation.pos not in used_qubits:
                    used_qubits.append(operation.pos)

            else:
                pos1 = min(operation.pos)
                pos2 = max(operation.pos)
                operationQlist[pos1].append(operation)
                for j in range(pos1 + 1, pos2 + 1):
                    operationQlist[j].append(None)
            
                for pos in operation.pos:
                    if pos not in used_qubits:
                        used_qubits.append(pos)

                maxlayer = max([len(operationQlist[j]) for j in range(pos1, pos2 + 1)])
                for j in range(pos1, pos2 + 1):
                    layerj = len(operationQlist[j])
                    pos = layerj - 1
                    if not layerj == maxlayer:
                        for i in range(abs(layerj - maxlayer)):
                            operationQlist[j].insert(pos, None)

        maxdepth = max([len(operationQlist[i]) for i in range(num)])

        for operations in operationQlist:
            operations.extend([None] * (maxdepth - len(operations)))

        # for m in self.measures.keys():
        #     if m not in used_qubits:
        #         used_qubits.append(m)
        used_qubits = np.sort(used_qubits)

        new_operationQlist = []
        for old_qi in range(len(operationQlist)):
            operations = operationQlist[old_qi]
            if old_qi in used_qubits:
                new_operationQlist.append(operations)

        lc = np.array(new_operationQlist)
        lc = np.vstack((used_qubits, lc.T)).T
        self.circuit = lc
        self._used_qubits = list(used_qubits)

    def draw(self, width : int=4, return_str : bool=False):
        """
        Draw layered circuit using ASCII, print in terminal.

        Args:
            width (int): The width of each operation.
            return_str: Whether return the circuit string.
        """
        self.layered_circuit()
        operationQlist = self.circuit
        num = operationQlist.shape[0]
        depth = operationQlist.shape[1] - 1
        printlist = np.array([[""] * depth for i in range(2 * num)], dtype="<U30")
        
        reduce_map = dict(zip(operationQlist[:, 0], range(num)))
        reduce_map_inv = dict(zip(range(num), operationQlist[:, 0]))
        for l in range(depth):
            layeroperations = operationQlist[:, l + 1]
            maxlen = 1 + width
            for i in range(num):
                operation = layeroperations[i]
                if  isinstance(operation, Operation) and len(operation.pos) == 1:
                    printlist[i * 2, l] = operation.symbol
                    maxlen = max(maxlen, len(operation.symbol) + width)
                elif isinstance(operation, Barrier):
                    pos = [i for i in operation.pos if i in reduce_map.keys()]
                    q1 = reduce_map[min(pos)]
                    q2 = reduce_map[max(pos)]
                    printlist[2 * q1:2 * q2 + 1, l] = "||"
                    maxlen = max(maxlen, len("||"))
                elif isinstance(operation, Operation) and len(operation.pos) > 1:
                    q1 = reduce_map[min(operation.pos)]
                    q2 = reduce_map[max(operation.pos)]
                    printlist[2 * q1 + 1:2 * q2, l] = "|"
                    printlist[q1 * 2, l] = "#"
                    printlist[q2 * 2, l] = "#"
                    if  hasattr(operation, "ctrls"): #Controlled-Multiqubit operation
                        for ctrl in operation.ctrls:
                            printlist[reduce_map[ctrl] * 2, l] = "*"
                        
                        if operation._targ_name == "SWAP":
                            printlist[reduce_map[operation.targs[0]] * 2, l] = "x"
                            printlist[reduce_map[operation.targs[1]] * 2, l] = "x"
                        else:
                            tq1 = reduce_map[min(operation.targs)]
                            tq2 = reduce_map[max(operation.targs)]
                            printlist[tq1 * 2, l] = "#"
                            printlist[tq2 * 2, l] = "#"
                            if tq1 + tq2 in [reduce_map[ctrl] * 2 for ctrl in operation.ctrls]:
                                printlist[tq1 + tq2, l] = "*" + operation.symbol
                            else:
                                printlist[tq1 + tq2, l] = operation.symbol
                            maxlen = max(maxlen, len(operation.symbol) + width)
                                
                    else: #Multiqubit operation
                        if operation.name == "SWAP":
                            printlist[q1 * 2, l] = "x"
                            printlist[q2 * 2, l] = "x"

                        else:
                            printlist[q1 + q2, l] = operation.symbol
                            maxlen = max(maxlen, len(operation.symbol) + width)

            printlist[-1, l] = maxlen

        circuitstr = []
        for j in range(2 * num - 1):
            if j % 2 == 0:
                linestr = ("q[%d]" % (reduce_map_inv[j // 2])).ljust(6) + "".join(
                    [printlist[j, l].center(int(printlist[-1, l]), "-") for l in range(depth)])
                if reduce_map_inv[j // 2] in self.measures.keys():
                    linestr += " M->c[%d]" % self.measures[reduce_map_inv[j // 2]]
                circuitstr.append(linestr)
            else:
                circuitstr.append("".ljust(6) + "".join(
                    [printlist[j, l].center(int(printlist[-1, l]), " ") for l in range(depth)]))
        circuitstr = "\n".join(circuitstr)
    
        if return_str:
            return circuitstr
        else:
            print(circuitstr+"\n")
   

    def to_qasm(self) -> str:
        """
        Convert the circuit to openqasm text.

        Returns: 
            openqasm text.
        """
        qasm = "OPENQASM 2.0;\ninclude \"qelib1.inc\";\n"
        qasm += "qreg q[%d];\n" % self.num
        qasm += "creg meas[%d];\n" % len(self.measures)
        for operation in self.operations:
            qasm += operation.to_qasm() + ";\n"

        for key in self.measures:
            qasm += "measure q[%d] -> meas[%d];\n" % (key, self.measures[key])

        self.openqasm = qasm
        return qasm


    def measure(self, pos: List[int], cbits: List[int] = []) -> None:
        """
        Measurement setting for experiment device.
        
        Args:
            pos: Qubits need measure.
            cbits: Classical bits keeping the measure results.
        """

        self.measures = dict(zip(pos, range(len(pos))))

        if cbits:
            if len(cbits) == len(self.measures):
                self.measures = dict(zip(pos, cbits))
            else:
                raise CircuitError("Number of measured bits should equal to the number of classical bits")


 