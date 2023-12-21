
import numpy as np
import matplotlib.pyplot as plt
from ..utils.basis import *
from ..exceptions import IqcsError
from .py_simulator import ptrace, permutebits

class SimuResult(object):
    def __init__(self, res_dict : dict):
        self._meta_data  =  res_dict
        
    def __getitem__(self, key:str):
        """
        Get meta_data of simulate results.
        Args:
            `"statevector"`: full state vector
            `"counts"`: sampled  bitstring counts
            `"pauli_expects"`: pauli expectations of input paulistrings
        """
        return self._meta_data[key]

        
    def plot_probabilities(self, full: bool=False, reverse_basis: bool=False, sort:bool=None):
        """
        Plot the probabilites of measured qubits
        """
        import matplotlib.pyplot as plt
        if self._meta_data["counts"]:
            counts = self._meta_data["counts"]
            total_counts = sum(counts.values())
            probabilities = {}
            for key in self._meta_data["counts"]:
                probabilities[key] = counts[key]/total_counts

            
            bitstrs = list(probabilities.keys())
            probs = list(probabilities.values())
            plt.figure()
            plt.bar(range(len(probs)), probs, tick_label = bitstrs)
            plt.xticks(rotation=70)
            plt.ylabel("probabilities")
            
        elif len(self._meta_data["statevector"]) > 0:
            psi =  self._meta_data["statevector"]
            num = int(np.log2(psi.shape[0]))
            measures = self._meta_data["measures"]
            psi = permutebits(psi, range(num)[::-1])
            probabilities = []
            if measures:
                probabilities = ptrace(psi, list(measures.keys()))
                probabilities = permutebits(probabilities, list(measures.values()))
                num = len(measures)
            else:
                probabilities = np.abs(psi)**2
            probs = probabilities
            inds = range(len(probs))
            if not full:
                inds = np.where(probabilities > 1e-14)[0]
                probs = probabilities[inds]

            basis=np.array([bin(i)[2:].zfill(num) for i in inds])
            if reverse_basis:
                basis=np.array([bin(i)[2:].zfill(num)[::-1] for i in inds])

            if sort == "ascend":
                orders = np.argsort(probs)
                probs = probs[orders]
                basis = basis[orders]
            elif sort == "descend":
                orders = np.argsort(probs)
                probs = probs[orders][::-1]
                basis = basis[orders][::-1]

            plt.figure()
            plt.bar(inds, probs, tick_label=basis)
            plt.xticks(rotation=70)
            plt.ylabel("probabilities")
        else:
            raise ValueError("No data for ploting")


