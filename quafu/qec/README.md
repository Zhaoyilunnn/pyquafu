# Design of the `qec` Module

## Overview

This document explains how the `qec` (Quantum Error Correction) module works in the `pyquafu` library. The main goal is to give users a way to build, test, and study quantum error correction codes. These codes help protect quantum information from noise. The module lets users try different codes and see how they perform when there is noise.

## Module Structure

The `qec` module has two main parts. One part is for quantum error correction codes. The other part is for decoders. There is also a file with base classes that set the rules for how codes and decoders should work.

The `quafu.qec.codes` part holds different quantum error correction codes. The first codes to be included are qLDPC codes and surface codes. The `quafu.qec.decoders` part has different ways to decode errors. Decoders use syndrome data to guess what errors happened. The `quafu.qec.base` file has base classes for codes and decoders. These base classes make sure that all codes and decoders use the same kind of interface. This makes it easier to add new codes or decoders later.

## Noise Model

A noise model is important when testing quantum error correction codes. The `qec` module will have a `noise_model.py` file. This file will let users add noise to their tests. The design of this noise model is based on the circuit-level noise model from [Stim](https://github.com/quantumlib/Stim/blob/main/doc/getting_started.ipynb). In this model, noise is added to all data qubits at certain points in the circuit. For example, a depolarizing channel can be added to all data qubits before each round of syndrome measurements, e.g., `before_round_data_depolarization`.

The noise channels use classes from `quafu.elements.noise`. Some of these classes are `Depolarizing`, `BitFlip`, and `Dephasing`. These classes help users add different types of noise to their tests.

## Open Questions and Future Work

Some problems still need to be solved. One problem is how to map qLDPC codes to physical layouts. The [`qLDPC`](https://github.com/oscarhiggott/PyMatching) library can help build qLDPC codes. But it is not clear how to map the check matrices ($C_X$ and $C_Z$) to real qubit positions on a 2D grid. This mapping is important for surface codes.

Another problem is how to make matching graphs for decoders. Decoders like [`PyMatching`](https://github.com/oscarhiggott/PyMatching) need a matching graph with edge weights. These weights should show how likely different errors are. It is not clear how to make this graph from a given circuit and the noise model. This is something that needs more work.
