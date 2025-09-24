"""

d = 3:
          X
    0---1---2
 Z  |   |   |
    3---4---5
    |   |   | Z
    6---7---8
      X

d = 5:
            X       X
      0---1---2---3---4
   Z  |   |   |   |   |
      5---6---7---8---9
      |   |   |   |   | Z
      10--11--12--13--14
   Z  |   |   |   |   |
      15--16--17--18--19
      |   |   |   |   | Z
      20--21--22--23--24
        X       X
"""

from quafu.qec.codes.surface_code import (
    PlaquetteType,
    _get_start_qubit_idx_weight_four,
    _get_start_qubit_idx_weight_two,
    _get_data_qubits_in_touch_order,
    Basis,
)


def test_get_touch_order():
    ############# weight-2 plaquettes ################
    # d=3
    assert _get_data_qubits_in_touch_order(
        0, 3, Basis.Z, PlaquetteType.weight_two_left
    ) == [None, None, 0, 3]
    assert _get_data_qubits_in_touch_order(
        5, 3, Basis.Z, PlaquetteType.weight_two_right
    ) == [5, 8, None, None]
    assert _get_data_qubits_in_touch_order(
        6, 3, Basis.X, PlaquetteType.weight_two_top
    ) == [6, 7, None, None]
    assert _get_data_qubits_in_touch_order(
        1, 3, Basis.X, PlaquetteType.weight_two_bottom
    ) == [None, None, 1, 2]

    # d=5
    assert _get_data_qubits_in_touch_order(
        0, 5, Basis.Z, PlaquetteType.weight_two_left
    ) == [None, None, 0, 5]
    assert _get_data_qubits_in_touch_order(
        10, 5, Basis.Z, PlaquetteType.weight_two_left
    ) == [None, None, 10, 15]
    assert _get_data_qubits_in_touch_order(
        9, 5, Basis.Z, PlaquetteType.weight_two_right
    ) == [9, 14, None, None]
    assert _get_data_qubits_in_touch_order(
        19, 5, Basis.Z, PlaquetteType.weight_two_right
    ) == [19, 24, None, None]
    assert _get_data_qubits_in_touch_order(
        20, 5, Basis.X, PlaquetteType.weight_two_top
    ) == [20, 21, None, None]
    assert _get_data_qubits_in_touch_order(
        22, 5, Basis.X, PlaquetteType.weight_two_top
    ) == [22, 23, None, None]
    assert _get_data_qubits_in_touch_order(
        1, 5, Basis.X, PlaquetteType.weight_two_bottom
    ) == [None, None, 1, 2]
    assert _get_data_qubits_in_touch_order(
        3, 5, Basis.X, PlaquetteType.weight_two_bottom
    ) == [None, None, 3, 4]


def test_get_start_qubit_idx():
    ################ weight-2 plaquettes ################
    # d=3
    assert _get_start_qubit_idx_weight_two(0, PlaquetteType.weight_two_left, 3) == 0
    assert _get_start_qubit_idx_weight_two(0, PlaquetteType.weight_two_right, 3) == 5
    assert _get_start_qubit_idx_weight_two(0, PlaquetteType.weight_two_top, 3) == 6
    assert _get_start_qubit_idx_weight_two(0, PlaquetteType.weight_two_bottom, 3) == 1

    # d=5
    assert _get_start_qubit_idx_weight_two(0, PlaquetteType.weight_two_left, 5) == 0
    assert _get_start_qubit_idx_weight_two(1, PlaquetteType.weight_two_left, 5) == 10
    assert _get_start_qubit_idx_weight_two(0, PlaquetteType.weight_two_right, 5) == 9
    assert _get_start_qubit_idx_weight_two(1, PlaquetteType.weight_two_right, 5) == 19
    assert _get_start_qubit_idx_weight_two(0, PlaquetteType.weight_two_top, 5) == 20
    assert _get_start_qubit_idx_weight_two(1, PlaquetteType.weight_two_top, 5) == 22
    assert _get_start_qubit_idx_weight_two(0, PlaquetteType.weight_two_bottom, 5) == 1
    assert _get_start_qubit_idx_weight_two(1, PlaquetteType.weight_two_bottom, 5) == 3

    ################ weight-4 plaquettes ################
    # d=3
    assert _get_start_qubit_idx_weight_four(0, 3, Basis.X) == 0
    assert _get_start_qubit_idx_weight_four(1, 3, Basis.X) == 4
    assert _get_start_qubit_idx_weight_four(0, 3, Basis.Z) == 1
    assert _get_start_qubit_idx_weight_four(1, 3, Basis.Z) == 3

    # d=5
    assert _get_start_qubit_idx_weight_four(0, 5, Basis.X) == 0
    assert _get_start_qubit_idx_weight_four(1, 5, Basis.X) == 2
    assert _get_start_qubit_idx_weight_four(2, 5, Basis.X) == 6
    assert _get_start_qubit_idx_weight_four(3, 5, Basis.X) == 8
    assert _get_start_qubit_idx_weight_four(4, 5, Basis.X) == 10
    assert _get_start_qubit_idx_weight_four(5, 5, Basis.X) == 12
    assert _get_start_qubit_idx_weight_four(6, 5, Basis.X) == 16
    assert _get_start_qubit_idx_weight_four(7, 5, Basis.X) == 18
    assert _get_start_qubit_idx_weight_four(0, 5, Basis.Z) == 1
    assert _get_start_qubit_idx_weight_four(1, 5, Basis.Z) == 3
    assert _get_start_qubit_idx_weight_four(2, 5, Basis.Z) == 5
    assert _get_start_qubit_idx_weight_four(3, 5, Basis.Z) == 7
    assert _get_start_qubit_idx_weight_four(4, 5, Basis.Z) == 11
    assert _get_start_qubit_idx_weight_four(5, 5, Basis.Z) == 13
    assert _get_start_qubit_idx_weight_four(6, 5, Basis.Z) == 15
    assert _get_start_qubit_idx_weight_four(7, 5, Basis.Z) == 17
