# TODO List for Modeling QEC Concepts in `quafu/qec/`

## 1. Implement Data Structures in `quafu/qec/`

- [ ] Create a `Location` class to represent (qubit, time_step) pairs.
- [ ] Create a `Detector` class to encapsulate measurement checks and detection logic.
- [ ] Create a `DetectingRegion` class to represent the set of locations and Pauli types a detector is sensitive to.
- [ ] Implement stabilizer flow logic (e.g., as a function or method) to propagate detecting regions through circuit operations.

## 2. Integrate with `QuantumCircuit`

- [ ] Add methods to `QuantumCircuit` to extract locations from the circuit.
- [ ] Add support for registering and managing detectors and their detecting regions.
- [ ] (Optional) Add methods to simulate or visualize stabilizer flow through the circuit.

## 3. Testing and Validation

- [ ] Write unit tests for each new class and method.
- [ ] Create example circuits and demonstrate how locations, detectors, and detecting regions are constructed and evolve.

## 4. Documentation

- [ ] Document the new classes and methods.
- [ ] Update the README in `quafu/qec/` to explain the new QEC modeling approach.

## 5. (Optional) Extend for Specific Codes

- [ ] Implement utilities for common QEC codes (e.g., surface code, qLDPC) to automatically generate detectors and detecting regions.
