from quafu.qec.location import get_all_locations_from_circ, Location
from quafu.circuits.quantum_circuit import QuantumCircuit


def test_get_all_locations_empty_circuit():
    qc = QuantumCircuit(2)
    locations = get_all_locations_from_circ(qc)
    # No instructions, so each qubit should have one location
    assert len(locations) == 2
    for loc in locations:
        assert loc.op_before is None
        assert loc.op_after is None


def test_get_all_locations_single_gate():
    qc = QuantumCircuit(2)
    qc.h(0)
    locations = get_all_locations_from_circ(qc)
    # Qubit 0: before H, after H
    # Qubit 1: no operations
    qubit0_locs = [loc for loc in locations if loc.qubit == 0]
    qubit1_locs = [loc for loc in locations if loc.qubit == 1]
    assert len(qubit0_locs) == 2
    assert len(qubit1_locs) == 1
    assert qubit1_locs[0].op_before is None
    assert qubit1_locs[0].op_after is None


def test_get_all_locations_multiple_gates():
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.x(0)
    qc.h(1)
    locations = get_all_locations_from_circ(qc)
    qubit0_locs = [loc for loc in locations if loc.qubit == 0]
    qubit1_locs = [loc for loc in locations if loc.qubit == 1]
    # Qubit 0: before H, between H and X, after X
    assert len(qubit0_locs) == 3
    # Qubit 1: before H, after H
    assert len(qubit1_locs) == 2
    assert qubit0_locs[0].op_before is None
    assert qubit0_locs[0].op_after is not None
    assert qubit0_locs[0].op_after.name == "H"
    assert qubit0_locs[1].op_before is not None
    assert qubit0_locs[1].op_before.name == "H"
    assert qubit0_locs[1].op_after is not None
    assert qubit0_locs[1].op_after.name == "X"
    assert qubit0_locs[2].op_before is not None
    assert qubit0_locs[2].op_before.name == "X"
    assert qubit0_locs[2].op_after is None


def test_get_all_locations_measure_and_gate():
    qc = QuantumCircuit(2)
    qc.h(0)
    qc.measure([0, 1], [0, 1])
    locations = get_all_locations_from_circ(qc)
    qubit0_locs = [loc for loc in locations if loc.qubit == 0]
    qubit1_locs = [loc for loc in locations if loc.qubit == 1]
    # Qubit 0: before H, between H and measure, after measure
    assert len(qubit0_locs) == 3
    # Qubit 1: before measure, after measure
    assert len(qubit1_locs) == 2
    assert qubit0_locs[0].op_before is None
    assert qubit0_locs[0].op_after is not None
    assert qubit0_locs[0].op_after.name == "H"
    assert qubit0_locs[1].op_before is not None
    assert qubit0_locs[1].op_before.name == "H"
    assert qubit0_locs[1].op_after is not None
    assert qubit0_locs[1].op_after.name == "measure"
    assert qubit0_locs[2].op_before is not None
    assert qubit0_locs[2].op_before.name == "measure"
    assert qubit0_locs[2].op_after is None
    assert qubit1_locs[0].op_before is None
    assert qubit1_locs[0].op_after is not None
    assert qubit1_locs[0].op_after.name == "measure"
    assert qubit1_locs[1].op_before is not None
    assert qubit1_locs[1].op_before.name == "measure"
    assert qubit1_locs[1].op_after is None
