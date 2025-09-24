"""
Implementation of the surface code

The modeling of plaquettes is inspired by tqec: https://github.com/tqec/tqec/
"""

import enum
from dataclasses import dataclass, field
from typing import Any, List
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


class Basis(enum.Enum):
    X = "X"
    Z = "Z"


def _get_data_qubits_in_touch_order(
    start_idx: int, d: int, basis: Basis, plaquette_type: PlaquetteType
) -> list[Any]:
    """
    Args:
        start_idx (int): The index of the top-left data qubit of the stabilizer.
            For weight-4 plaquettes this is the top-left qubit of the 2x2 block.
            For weight-2 plaquettes, this is the index of the left or top qubit.
        d (int): The distance of the surface code.
        basis (Basis): Basis.X or Basis.Z, indicating the type of stabilizer.
        plaquatte_type (PlaquetteType): The type of the plaquette.

    Returns:
        list[int]: The indices of the data qubits in the order they are touched by the
                   ancilla qubit.
    """
    if plaquette_type == PlaquetteType.weight_four:
        if basis == Basis.X:
            return [
                start_idx,
                start_idx + 1,
                start_idx + d,
                start_idx + d + 1,
            ]
        elif basis == Basis.Z:
            return [
                start_idx,
                start_idx + d,
                start_idx + 1,
                start_idx + d + 1,
            ]
        else:
            raise ValueError("basis must be Basis.X or Basis.Z")
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


def _get_start_qubit_idx_weight_four(plaquette_idx: int, d: int, basis: Basis) -> int:
    """
    Get the index of the top-left data qubit of a weight-4 plaquette.
    """
    assert 0 <= plaquette_idx < (d - 1) ** 2 // 2
    plaquettes_per_row = (d - 1) // 2
    row_id = plaquette_idx // plaquettes_per_row
    col_id = plaquette_idx % plaquettes_per_row
    parity = row_id % 2
    if basis == Basis.X:
        return (d * row_id) + col_id * 2 + parity
    elif basis == Basis.Z:
        return (d * row_id) + col_id * 2 + (1 - parity)


@dataclass
class Plaquette:
    """
    Type of stabilizer plaquette.
    Simply consists of an ancilla qubit id, plaquette_type, plaquette_idx, and basis.
    """

    ancilla_idx: int
    plaquette_type: PlaquetteType
    plaquette_idx: int
    basis: Basis  # X or Z
    d: int
    touch_order: List[Any] = field(init=False, default_factory=list)

    def __post_init__(self):
        """
        Get the indices of the data qubits in the order they are touched by the
        ancilla qubit.
        """
        if self.plaquette_type == PlaquetteType.weight_four:
            start_idx = _get_start_qubit_idx_weight_four(
                self.plaquette_idx, self.d, self.basis
            )
        else:
            start_idx = _get_start_qubit_idx_weight_two(
                self.plaquette_idx, self.plaquette_type, self.d
            )
        self.touch_order = _get_data_qubits_in_touch_order(
            start_idx, self.d, self.basis, self.plaquette_type
        )
        assert len(self.touch_order) == 4


class SurfaceCode:
    def __init__(self, d: int):
        """
        Initialize the SurfaceCode object.

        Args:
            d (int): The distance of the surface code.
        """
        self.d = d
        self.num_data_qubits = d**2
        self.num_ancilla_qubits = d**2 - 1
        self.num_qubits = self.num_data_qubits + self.num_ancilla_qubits
        self.plaquettes: List[Plaquette] = []
        self._create_plaquettes()

    def _create_plaquettes(self):
        """
        Create all plaquettes and store them in self.plaquettes.
        """
        num_weight_four_plaquettes = (self.d - 1) ** 2 // 2
        num_weight_two_plaquettes = (self.d - 1) // 2

        ancilla_idx = 0
        plaquette_metadata = [
            (["weight_four"], Basis.X, num_weight_four_plaquettes),
            (["weight_four"], Basis.Z, num_weight_four_plaquettes),
            (
                ["weight_two_left", "weight_two_right"],
                Basis.Z,
                num_weight_two_plaquettes,
            ),
            (
                ["weight_two_top", "weight_two_bottom"],
                Basis.X,
                num_weight_two_plaquettes,
            ),
        ]
        for types, pauli, num_plaquettes in plaquette_metadata:
            for plaquette_type in [getattr(PlaquetteType, t) for t in types]:
                for i in range(num_plaquettes):
                    plaquette = Plaquette(
                        ancilla_idx=ancilla_idx,
                        plaquette_type=plaquette_type,
                        plaquette_idx=i,
                        basis=pauli,
                        d=self.d,
                    )
                    self.plaquettes.append(plaquette)
                    ancilla_idx += 1

        assert ancilla_idx == self.num_ancilla_qubits

    def create_circuit(self) -> QuantumCircuit:
        """
        Create the quantum circuit for the surface code.

        Returns:
            QuantumCircuit: The constructed quantum circuit.
        """
        ancilla_idxes = [plaquette.ancilla_idx for plaquette in self.plaquettes]
        data_idxes = list(range(self.num_data_qubits))

        qc = QuantumCircuit(self.num_qubits, self.num_qubits)

        # 1. Initialization
        # All data qubits to |0⟩
        qc.reset(data_idxes)

        # 2. d rounds of syndrome measurements
        for _ in range(self.d):
            # Reset ancilla qubits
            qc.reset(ancilla_idxes)
            for plaquette in self.plaquettes:
                if plaquette.basis == Basis.X:
                    qc.h(plaquette.ancilla_idx)
            # CNOTs
            for i in range(4):
                for plaquette in self.plaquettes:
                    if plaquette.basis == Basis.X:
                        if plaquette.touch_order[i] is not None:
                            target = plaquette.touch_order[i]
                            qc.cx(plaquette.ancilla_idx, target)
                    elif plaquette.basis == Basis.Z:
                        if plaquette.touch_order[i] is not None:
                            control = plaquette.touch_order[i]
                            qc.cx(control, plaquette.ancilla_idx)
            # Measurement
            for plaquette in self.plaquettes:
                if plaquette.basis == Basis.X:
                    qc.h(plaquette.ancilla_idx)
            qc.measure(ancilla_idxes, ancilla_idxes)

        # 3. Measure all data qubits in Z basis
        qc.measure(data_idxes, data_idxes)

        return qc
