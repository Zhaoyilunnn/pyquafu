# Design of the `qec` Module

## Overview

This document outlines the design of the `qec` (Quantum Error Correction) module for the `pyquafu` library. The primary goal of this module is to provide a flexible and extensible framework for building, simulating, and analyzing the performance of quantum error correction codes under various noise models.

## Module Structure

The `qec` module is organized into two main submodules, `codes` and `decoders`, with foundational abstract classes defined in `base.py`.

* **`quafu.qec.codes`**: This submodule will contain implementations of various quantum error correction codes. The design is intended to be general, starting with support for qLDPC codes and surface codes.
* **`quafu.qec.decoders`**: This submodule will house different decoding algorithms. These decoders will process syndrome information to infer the most likely errors that have occurred.
* **`quafu.qec.base`**: This file defines the abstract base classes for codes and decoders, ensuring a consistent API and promoting extensibility for future contributions.

## Noise Model

A crucial component for evaluating QEC codes is a realistic noise model. We will implement a `noise_model.py` module within the `qec` directory.

The design of this noise model is inspired by Stim's circuit-level noise model (see [Stim's getting started guide](https://github.com/quantumlib/Stim/blob/main/doc/getting_started.ipynb)). The core idea is to apply noise channels to all data qubits at specific stages of the quantum circuit's execution. For example, a `before_round_data_depolarization` channel could be applied to all data qubits before each round of syndrome measurements.

To implement the noise channels, we will reuse the existing classes available in `quafu.elements.noise`, such as `Depolarizing`, `BitFlip`, and `Dephasing`.

## Open Questions and Future Work

There are several areas that require further research and development:

1. **Mapping qLDPC Codes to Physical Layouts**: We plan to leverage the [`qLDPC`](https://github.com/oscarhiggott/PyMatching) library for the construction of general qLDPC codes. A key challenge is to establish a clear and systematic mapping from the abstract check matrices ($C_X$ and $C_Z$) of a given code to the physical qubit coordinates on a 2D grid, which is particularly important for implementing surface codes.

2. **Generating Matching Graphs for Decoders**: For integration with powerful decoders like [`PyMatching`](https://github.com/oscarhiggott/PyMatching), it is necessary to construct a matching graph with appropriate edge weights that reflect the error probabilities. It is currently an open question how to systematically generate this weighted graph from a given circuit structure and our defined circuit-level noise model.

