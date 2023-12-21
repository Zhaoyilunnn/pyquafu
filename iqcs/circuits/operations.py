
from abc import ABC, abstractmethod

class Operation(ABC):
    @property
    @abstractmethod
    def name(self):
        """Unique string identifier for operation type."""
        raise NotImplementedError


    @property
    @abstractmethod
    def pos(self):
        """Qubits the operation acts on"""
        raise NotImplementedError
    

class Barrier(Operation):
    name = "barrier"
    def __init__(self, pos:list):
        self.name = "barrier"
        self.__pos = pos
        self.symbol = "||"

    @property
    def pos(self):
        return self.__pos

    @pos.setter
    def pos(self, pos):
        self.__pos = pos

    def __repr__(self):
        return f"{self.__class__.__name__}"

    def to_qasm(self):
        return "barrier " + ",".join(["q[%d]" % p for p in range(min(self.pos), max(self.pos)+1)])

class Delay(Operation):
    name = "delay"
    def __init__(self, pos : int, duration : int, unit="ns"):
        if isinstance(duration, int):
            self.duration = duration
        else:
            raise TypeError("duration must be int")
        self.unit=unit
        self.pos=[pos]
        self.symbol = "Delay(%d%s)" %(duration, unit)

    def __repr__(self):
        return f"{self.__class__.__name__}"

    def to_qasm(self):
        return "delay(%d%s) q[%d]" % (self.duration, self.unit, self.pos[0])

class Measure(object):
    name = "measure"
    def __init__(self, bitmap : dict):
        self.qbits = bitmap.keys()
        self.cbits = bitmap.values()