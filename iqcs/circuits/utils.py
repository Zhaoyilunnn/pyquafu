
import numpy as np
from typing import List, Union
from .parameters import Parameter, ParameterExpression

def reorder_matrix(matrix : np.ndarray, pos : List):
    """Reorder the input sorted matrix to the pos order """
    qnum = len(pos)
    dim = 2**qnum
    inds = np.argsort(pos)
    inds = np.concatenate([inds, inds+qnum])
    tensorm = np.reshape(matrix, [2]*2*qnum)
    return np.transpose(tensorm, inds).reshape([dim, dim])

def extract_float(paras):
    paras_f = []
    for para in paras:
        if isinstance(para, float) or isinstance(para, int):
            paras_f.append(para)
        elif isinstance(para, Parameter) or isinstance(para, ParameterExpression):
            paras_f.append(para.get_value())
    return paras_f