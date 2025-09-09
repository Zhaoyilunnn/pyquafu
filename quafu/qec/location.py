from dataclasses import dataclass
from typing import Optional
from quafu.elements.instruction import Instruction
from quafu import QuantumCircuit


@dataclass
class Location:
    """
    Represents a location in a quantum circuit: a specific qubit between specific operations.

    References:
        - http://arxiv.org/abs/2302.02192
    """

    qubit: int
    op_before: Optional[Instruction] = None
    op_after: Optional[Instruction] = None


def get_all_locations_from_circ(circuit: QuantumCircuit) -> list[Location]:
    """
    Get all possible locations in a quantum circuit.

    Args:
        circuit (QuantumCircuit): The quantum circuit to analyze.

    Returns:
        list[Location]: A list of Location objects representing all possible locations in the circuit.
    """
    locations = []
    for qubit in range(circuit.num):
        ops_on_qubit = [op for op in circuit.instructions if qubit in op.pos]
        if not ops_on_qubit:
            # If no operations on this qubit, only one location exists
            locations.append(Location(qubit))
        else:
            # Location before the first operation
            locations.append(Location(qubit, op_before=None, op_after=ops_on_qubit[0]))
            # Locations between operations
            for i in range(len(ops_on_qubit) - 1):
                locations.append(
                    Location(
                        qubit, op_before=ops_on_qubit[i], op_after=ops_on_qubit[i + 1]
                    )
                )
            # Location after the last operation
            locations.append(Location(qubit, op_before=ops_on_qubit[-1], op_after=None))
    return locations
