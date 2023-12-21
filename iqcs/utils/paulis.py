# from ..exceptions import IqcsError
import scipy.sparse as sp
import numpy as np

IMat = sp.coo_matrix(np.array([[1., 0.],
                               [0., 1.]],dtype=complex))

XMat = sp.coo_matrix(np.array([[0., 1.],
                               [1., 0.]], dtype=complex))

YMat = sp.coo_matrix(np.array([[0., -1.j],
                               [1.j, 0.]], dtype=complex))

ZMat = sp.coo_matrix(np.array([[1., 0.],    
                               [0., -1.]], dtype=complex))

PauliMats = {"X": XMat, "Y":YMat, "Z":ZMat, "I":IMat}

class PauliOp(object):
    def __init__(self, paulis:str, coeff = 1.):
        paulist = paulis.split(" ")
        self.paulistr = ''
        self.pos = []
        for p in paulist:
            assert p[0] in "XYZ"
            self.paulistr += p[0]
            self.pos.append(int(p[1:]))
        self.coeff = coeff

    def __repr__(self):
        repstr = ""
        if self.coeff != 1.:
            repstr = str(self.coeff) + "*"
        for i in range(len(self.pos)):
            repstr += self.paulistr[i]
            repstr += str(self.pos[i])
            repstr += "*"
        return  repstr[:-1]
    
    def __str__(self):
        return self.__repr__()

    def __mul__(self, obj):
        pass

    def __rmul__(self, obj):
        pass
    
    def commutator(self, obj):
        pass
    
    def get_matrix(self, qnum, big_endian=False):
        if (qnum - 1 < max(self.pos)):
            raise ValueError("The support of the paulis exceed the total qubit number")

        pos = np.array(self.pos)
        if not big_endian:
            pos = qnum - 1 - pos
        inds = np.argsort(pos)
        iq = 0
        ip = 0
        mat = 1.
        while iq < qnum:
            if ip < len(pos):
                if iq == pos[inds[ip]]:
                    opstr = self.paulistr[inds[ip]]
                    mat = sp.kron(mat, PauliMats[opstr])
                    iq += 1
                    ip += 1
                else:
                    mat = sp.kron(mat,PauliMats["I"])
                    iq += 1
            else:
                mat = sp.kron(mat,PauliMats["I"])
                iq += 1

        return self.coeff * mat

class Hamiltonian(object):
    def __init__(self, paulis:list[PauliOp]):
        self.paulis = paulis
        self.matrix = None

    def __repr__(self):
        return "+".join([str(pauli) for pauli in self.paulis])

    def __str__(self):
        return self.__repr__()

    def get_matrix(self, qnum, big_endian=False):
        mat = 0.
        for pauli in self.paulis:
            mat += pauli.get_matrix(qnum, big_endian)   
        self.matrix = mat
        return mat
    

def intersec(a, b):
    inter = []
    aind = []
    bind = []
    for i in range(len(a)):
        for j in range(len(b)):
            if a[i] == b[j]:
                inter.append(a[i])
                aind.append(i)
                bind.append(j)
    
    return inter, aind, bind

def diff(a, b):
    diff = []
    aind = []
    for i in range(len(a)):
        if a[i] not in b:
            diff.append(a[i])
            aind.append(i)
    
    return diff, aind

def merge_paulis(obslist):
    measure_basis = []
    targ_basis = []
    for obs in obslist:
        if len(measure_basis) == 0:
            measure_basis.append(obs)
            targ_basis.append(len(measure_basis)-1)
        else:
            added = 0
            for mi in range(len(measure_basis)):
                measure_base = measure_basis[mi]
                interset, intobsi, intbasei = intersec(obs.pos, measure_base.pos) 
                diffset, diffobsi = diff(obs.pos, measure_base.pos)
                if not len(interset) == 0:
                    if all(np.array(list(obs.paulistr))[intobsi] == np.array(list(measure_base.paulistr))[intbasei]):
                        measure_base.paulistr += "".join(np.array(list(obs.paulistr))[diffobsi])
                        measure_base.pos.extend(diffset)
                        targ_basis.append(mi)
                        added = 1
                        break
                else:
                    measure_base.paulistr += obs.paulistr
                    measure_base.pos.extend(obs.pos)
                    targ_basis.append(mi)
                    added = 1
                    break

            if not added: 
                measure_basis.append(obs)
                targ_basis.append(len(measure_basis)-1)

    return measure_basis, targ_basis



if __name__ == "__main__":
    op1 = PauliOp("X2 Z0")
    mat1 = op1.get_matrix(4).toarray()
    mat2 = sp.kron(sp.kron(PauliMats["Z"], PauliMats["I"]), PauliMats["X"])
    mat2 = sp.kron(mat2, PauliMats["I"])
    mat2 = mat2.toarray()
    print(np.linalg.norm(mat1-mat2))