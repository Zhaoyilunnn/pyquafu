from quafu.qec.codes.surface_code import (
    PlaquetteType,
    _get_start_qubit_idx_weight_four,
    _get_start_qubit_idx_weight_two,
)


def test_get_start_qubit_idx():
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
    assert _get_start_qubit_idx_weight_four(0, 3, "X") == 0
    assert _get_start_qubit_idx_weight_four(1, 3, "X") == 4
    assert _get_start_qubit_idx_weight_four(0, 3, "Z") == 1
    assert _get_start_qubit_idx_weight_four(1, 3, "Z") == 3

    # d=5
    assert _get_start_qubit_idx_weight_four(0, 5, "X") == 0
    assert _get_start_qubit_idx_weight_four(1, 5, "X") == 2
    assert _get_start_qubit_idx_weight_four(2, 5, "X") == 6
    assert _get_start_qubit_idx_weight_four(3, 5, "X") == 8
    assert _get_start_qubit_idx_weight_four(4, 5, "X") == 10
    assert _get_start_qubit_idx_weight_four(5, 5, "X") == 12
    assert _get_start_qubit_idx_weight_four(6, 5, "X") == 16
    assert _get_start_qubit_idx_weight_four(7, 5, "X") == 18
    assert _get_start_qubit_idx_weight_four(0, 5, "Z") == 1
    assert _get_start_qubit_idx_weight_four(1, 5, "Z") == 3
    assert _get_start_qubit_idx_weight_four(2, 5, "Z") == 5
    assert _get_start_qubit_idx_weight_four(3, 5, "Z") == 7
    assert _get_start_qubit_idx_weight_four(4, 5, "Z") == 11
    assert _get_start_qubit_idx_weight_four(5, 5, "Z") == 13
    assert _get_start_qubit_idx_weight_four(6, 5, "Z") == 15
    assert _get_start_qubit_idx_weight_four(7, 5, "Z") == 17
