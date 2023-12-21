from iqcs.circuits.quantum_circuit import QuantumCircuit

from typing import List
from iqcs.exceptions import CircuitError, ServerError
import requests
import json
from iqcs.utils.convertors import circuit_to_json
from typing import List
from .client import *
from collections import OrderedDict
import numpy as np
from .utils.basis import *
from .exceptions import IqcsError
from copy import deepcopy
from .account import load_token
from .utils.paulis import PauliOp, merge_paulis
from iqcs.circuits.gatelib import RY, RX
from iqcs.circuits.operations import Barrier

class QuantumExperiment(object):
    """
    Quantum experiment
    """
    def __init__(self, qc: QuantumCircuit, name="circuit0",shots=1024, compile_circuit=False, group="", backend="simulation", measure_axis=[]):
        self._qc = qc
        self._origates = qc.operations
        self.name = name
        self.shots = shots
        self.compile_circuit = compile_circuit
        self.group = group
        self.backend = backend
        self.measure_axis = measure_axis

    @property
    def qc(self):
        if self.measure_axis:
            self._qc.operations = deepcopy(self._origates)
            self._qc << (Barrier(self.measure_axis.pos))
            for ax, pos in zip(self.measure_axis.paulistr, self.measure_axis.pos):
                if ax == "X":
                    self._qc << (RY(pos, -np.pi / 2))
                elif ax== "Y":
                    self._qc << (RX(pos, np.pi / 2))
        return self._qc
    
    @qc.setter
    def qc(self, qc):
        self._qc = qc

    def config(self, **kwargs):
        for attr, value in kwargs.items():
            if hasattr(self, attr):
                setattr(self, attr, value)
            else:
                raise ValueError("No such attribute")
    
    def draw_circuit(self):
        self.qc.draw()

class ExperimentResult(object):
    def __init__(self, job_id, name='', counts={}, status='', cqmap={}):
        self.job_id = job_id
        self.status = status
        self.cqmap = cqmap
        self.name = name
        self.counts = {}
        if counts:
            self.init_with_counts(counts)

    def init_with_counts(self, counts):
        self.counts = counts
        self._logical_counts = counts
    
        if self.cqmap:
            self._logical_counts = OrderedDict(sorted(self.counts.items(), key=lambda s: s[0]))
            cbits = list(self.cqmap.values())
            for key, values in self.counts.items():
                newkey = "".join([key[i] for i in cbits])
                self.logicalq_res[newkey] = values

        total_counts = sum(self.counts.values())
        self.probabilities = {} 
        for key in self.counts:
            self.probabilities[key] = self.counts[key]/total_counts

    def retrieve(self):
        data = {'job_id':self.job_id}
        r = requests.post(query_url ,data=data)
        res = json.loads(r.content)
        if res['code'] == "0":
            # print("status: ", res['data']['status'])
            self.status = res['data']['status']
            if res['data']['status'] == "DONE":
                xy_data = res["data"]["result"]["data"]
                raw_data = { p['X'] : p['Y'] for p in xy_data} 
                self.init_with_counts(raw_data)
        else:
            raise ServerError(res["errMsg"]+res["message"])
    
    def expectZ(self, pos):
        return measure_obs(pos, self._logical_counts)
    
    def plot_probabilities(self):
        import matplotlib.pyplot as plt
        bitstrs = list(self.probabilities.keys())
        probs = list(self.probabilities.values())
        plt.figure()
        plt.bar(range(len(probs)), probs, tick_label = bitstrs)
        plt.xticks(rotation=70)
        plt.ylabel("probabilities")

class ResultGroup(object):
    def __init__(self, res_list: List[ExperimentResult]=None, name='', merge_info:List=None):
        self.res_list = []
        self.merge_info = []
        if not res_list == None:
            self.res_list = res_list
        if not merge_info == None:
            self.merge_info = merge_info

        if self.merge_info:
            self.obs_list = self.merge_info[0]
        self.name = name
        self.obs_expct = []

    def check_status(self, verbose=False):
        done = True
        self.retrieve()
        if verbose:
            print((" "*5).join(["job_id".ljust(16), "job_name".ljust(10),   "status".ljust(10)]))
        for res in self.res_list:
            if verbose:
                name = res.name if res.name else "Untitled"
                print((" "*5).join([(res.job_id[:10]+"...").ljust(16), ("%s" %name).ljust(10), ("%s" %res.status).ljust(10)]))
            if not res.status == "DONE":
                done = False
        return done
    
    def retrieve(self):
        for res in self.res_list:
            res.retrieve()
        
    def expectZ(self, pos):
        if self.check_status():
            expcts = []
            for res in self.res_list:
                expcts.append(res.expectZ(pos))
            return expcts
        else:
            raise IqcsError("The job is not completed")
        
    def expectZ_mapped(self, pos_map):
        if self.check_status():
            expcts = []
            pos, rmap = pos_map
            for i in range(len(pos)):
                expcts.append(self.res_list[rmap[i]].expectZ(pos[i]))
            return expcts
        else:
            raise IqcsError("The job is not completed")
    
    def calculate_obs(self):
        #TODO:add coefficients
        if len(self.obs_list) == 0:
            raise IqcsError("Corresponding job has no observerable measured")
        pos_map = []
        pos_map.append([self.obs_list[i].pos for i in range(len(self.obs_list))])
        pos_map.append(self.merge_info[1])
        expcts = self.expectZ_mapped(pos_map)
        self.obs_expct = expcts

    def add_result(self, res):
        self.res_list.append(res)

def check_validity(exp: QuantumExperiment):
    num_backend = 0 
    if exp.backend == "computer":
        num_backend = 10
    elif exp.backend == "simulation":
        num_backend = 30
    if exp.qc.num > num_backend:
            raise CircuitError("The qubit number %d is too large for backend %s which has %d qubits" %(exp.qc.num, exp.backend, num_backend))

 
def submit_job(qexp: QuantumExperiment, wait:bool=True, obslist:List[PauliOp]=[]):
    token = load_token()
    measures = list(qexp.qc.measures.keys())
    if len(obslist) == 0: #single experiment
        data = {'experiments':[{'instructions':circuit_to_json(qexp.qc)}], 'schema_version': '1.1.0', 'type': qexp.backend, 'qobj_id':'8d7c9805-8764-48a7-b0ac-7c81a94459c8','device_id': '5ae875670f020500393162ad','token':token, 'config': {
            'shots': qexp.shots,
            'n_qubits': qexp.qc.num,
            'memory_slots':len(qexp.qc.measures), 
            'wait': wait
            }}

        data = json.dumps(data)
        #check_validity
        r = requests.post(validate_url, data=data)
        check_res = json.loads(r.content)
        # print(str(check_res))
        if check_res['code'] == '0':
            #execute
            r = requests.post(submit_url, data=data)
            res = json.loads(r.content)
            if wait:
                if res['code'] == "0":
                    xy_data = res["data"]["result"]["data"]
                    raw_data = { p['X'] : p['Y'] for p in xy_data}
                    return  ExperimentResult(res['data']['job_id'], counts = raw_data, status=res['data']['status'])
                else:
                    raise ServerError(res["errMsg"]+res["message"])
            else:
                if res['code'] == '0':
                    return ExperimentResult(res['data']['job_id'])
                
                else:
                    raise ServerError(res['errMsg']+res["message"])
        else:
            raise CircuitError(check_res["errMsg"]+check_res["message"])

    else: #group experiment for measure observables
        for obs in obslist:
            for p in obs.pos:
                if p not in measures:
                    raise CircuitError("Qubit %d in observer %s is not measured." % (p, obs))
   
        measure_axes, targlist = merge_paulis(obslist)

        print("Job start, need measured in ", measure_axes)
        group_res = ResultGroup(merge_info=[obslist, targlist])
        for measure_axis in measure_axes:
            qexp.config(measure_axis=measure_axis)
            res = submit_job(qexp, wait)
            group_res.add_result(res)
        if wait:
            group_res.calculate_obs()
        return group_res

def submit_jobs(exp_list: List[QuantumExperiment], group_name='', wait=True):
    group_res = ResultGroup(name=group_name)
    for exp in exp_list:
        group_res.add_result(submit_job(exp, wait=wait))

    return group_res


