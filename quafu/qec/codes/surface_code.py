"""
Implementation of the surface code
"""

import enum
from dataclasses import dataclass
from typing import Any
from quafu.circuits import QuantumCircuit


class PlaquetteType(enum.Enum):
    """
    Currently only assume the top/bottom-boundary is X-type
    and left/right-boundary is Z-type.
    """

    weight_four = 1
    weight_two_top = 2
    weight_two_bottom = 3
    weight_two_left = 4
    weight_two_right = 5


def _get_data_qubits_in_touch_order(
    start_idx: int, d: int, basis: str, plaquette_type: PlaquetteType
) -> list[Any]:
    """
    Args:
        start_idx (int): The index of the top-left data qubit of the stabilizer.
            For weight-4 plaquettes this is the top-left qubit of the 2x2 block.
            For weight-2 plaquettes, this is the index of the left or top qubit.
        d (int): The distance of the surface code.
        basis (str): 'X' or 'Z', indicating the type of stabilizer.
        plaquatte_type (PlaquetteType): The type of the plaquette.

    Returns:
        list[int]: The indices of the data qubits in the order they are touched by the
                   ancilla qubit.
    """
    if plaquette_type == PlaquetteType.weight_four:
        if basis == "X":
            return [
                start_idx,
                start_idx + 1,
                start_idx + d,
                start_idx + d + 1,
            ]
        elif basis == "Z":
            return [
                start_idx,
                start_idx + d,
                start_idx + 1,
                start_idx + d + 1,
            ]
        else:
            raise ValueError("basis must be 'X' or 'Z'")
    # Else: weight-2 plaquettes
    # Now we assume the top/bottom-boundary is X-type
    # and left/right-boundary is Z-type.
    else:
        if plaquette_type == PlaquetteType.weight_two_top:
            return [start_idx, start_idx + 1, None, None]
        elif plaquette_type == PlaquetteType.weight_two_bottom:
            return [None, None, start_idx, start_idx + 1]
        elif plaquette_type == PlaquetteType.weight_two_left:
            return [None, None, start_idx, start_idx + d]
        elif plaquette_type == PlaquetteType.weight_two_right:
            return [start_idx, start_idx + d, None, None]


def _get_start_qubit_idx_weight_two(
    plaquette_idx: int, plaquette_type: PlaquetteType, d: int
) -> int:
    """
    Get the index of the top or left data qubit of a weight-2 plaquette.
    """
    assert 0 <= plaquette_idx < (d - 1) // 2
    assert plaquette_type is not PlaquetteType.weight_four
    if plaquette_type == PlaquetteType.weight_two_left:
        return plaquette_idx * 2 * d
    elif plaquette_type == PlaquetteType.weight_two_right:
        return plaquette_idx * 2 * d + (2 * d - 1)
    elif plaquette_type == PlaquetteType.weight_two_top:
        return d * (d - 1) + plaquette_idx * 2
    elif plaquette_type == PlaquetteType.weight_two_bottom:
        return 1 + plaquette_idx * 2


def _get_start_qubit_idx_weight_four(plaquette_idx: int, d: int, basis: str) -> int:
    """
    Get the index of the top-left data qubit of a weight-4 plaquette.
    """
    assert 0 <= plaquette_idx < (d - 1) ** 2 // 2
    plaquettes_per_row = (d - 1) // 2
    row_id = plaquette_idx // plaquettes_per_row
    col_id = plaquette_idx % plaquettes_per_row
    parity = row_id % 2
    if basis == "X":
        return (d * row_id) + col_id * 2 + parity
    elif basis == "Z":
        return (d * row_id) + col_id * 2 + (1 - parity)
    else:
        raise ValueError("basis must be 'X' or 'Z'")


@dataclass
class Plaquette:
    """
    Type of stabilizer plaquette.
    Simply consists of an ancilla qubit id, plaquette_type, plaquette_idx, and basis.
    """

    ancilla_idx: int
    plaquette_type: PlaquetteType
    plaquette_idx: int
    basis: str  # 'X' or 'Z'


def canonical_surface_code_circuit(d: int):
    """
    This function is just a temporary implementation for testing purposes.
    It only supports single-logical-qubit rotated surface code, with `d` rounds of
    stabilizer measurements.

    Order of touching data qubits:
    - X stabilizers: Z-type
    - Z stabilizers: reversed N-type (top-left, bottom-left, top-right, bottom-right)

    Args:
        d (int): The distance of the surface code.
    """
    num_data_qubits = d**2
    num_ancilla_qubits = d**2 - 1
    num_qubits = num_data_qubits + num_ancilla_qubits

    # Number of weight-4 plaquettes per type (X or Z)
    num_weight_four_plaquettes = (d - 1) ** 2 // 2

    # Number of weight-2 plaquettes per type (top/bottom or left/right)
    num_weight_two_plaquettes = (d - 1) // 2

    ancilla_idx = 0

    plaquettes = []

    # Create plaquettes
    plaquette_types = [
        (["weight_four"], "X", num_weight_four_plaquettes),
        (["weight_four"], "Z", num_weight_four_plaquettes),
        (["weight_two_left", "weight_two_right"], "Z", num_weight_two_plaquettes),
        (["weight_two_top", "weight_two_bottom"], "X", num_weight_two_plaquettes),
    ]
    for types, pauli, num_plaquettes in plaquette_types:
        for plaquette_type in [getattr(PlaquetteType, t) for t in types]:
            for i in range(num_plaquettes):
                plaquette = Plaquette(ancilla_idx, plaquette_type, i, pauli)
                plaquettes.append(plaquette)
                ancilla_idx += 1

    assert ancilla_idx == num_ancilla_qubits

    # Create circuit
    qc = QuantumCircuit(num_qubits, num_qubits)

    # TODO: Add operations
