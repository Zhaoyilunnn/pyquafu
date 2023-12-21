#pragma once
#include "dagcircuit.hpp"

Circuit reduce_by_measure(Circuit const& circuit){
    auto dag = DAG::DAGCircuit(circuit);
    auto mind = dag.add_measures();
    dag.reverse();
    auto subdag = dag.find_connected({mind});
    subdag.reverse();
    auto gates = subdag.topological_sort();
    Circuit newcircuit(circuit.qubit_num());
    for (auto op : gates){
        if (op.name() != "measure"){
            newcircuit.add_gate(op);
        }
    }
    return newcircuit;
}

Circuit reduce_absorb(Circuit const& circuit){
    auto dag = DAG::DAGCircuit(circuit);
    dag.absorb_gates();
    auto newcircuit = dag.to_circuit();
    // std::cout << "gate num " << newcircuit.gates().size() << std::endl;
    return newcircuit ; 
}

Circuit reduce_merge(Circuit const& circuit, int max_block_size=5)
{
    auto dag = DAG::DAGCircuit(circuit);
    dag.merge_gates(max_block_size);
    auto newcircuit = dag.to_circuit();
    // std::cout << "gate num " << newcircuit.gates().size() << std::endl;
    return newcircuit ; 
}
