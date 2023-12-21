#pragma once
#include "operators.hpp"
#include "qasm.hpp"
#include "util.h"
#include <algorithm>


class Circuit{
    private:
        uint qubit_num_;    
        vector<QuantumOperator> gates_;
        uint max_targe_num_;
        std::map<int, int> measures_;
    public:
    Circuit();
    explicit Circuit(uint qubit_num);
    Circuit(uint qubit_num, vector<QuantumOperator> const &gates,  std::map<int, int> const& measures);
    Circuit(py::object const&pycircuit, bool get_full_mat=false, bool reverse=true); 

    void add_gate(QuantumOperator const &gate);
    uint qubit_num() const { return qubit_num_; }
    uint max_targe_num() const {return max_targe_num_;}
    vector<QuantumOperator>gates() const { return gates_; }
    std::map<int, int> measures() const {return measures_;}
    void set_measures(std::map<int, int> const measures) { measures_ = measures; }
    void set_gates(std::vector<QuantumOperator> const& gates){ gates_ = gates;}

    void print(){
        for (auto op : gates_){
            std::cout << op.name() << " ";
            for (auto p : op.positions()){
                std::cout << p << " ";
            }
            std::cout << std::endl;
        }
    }
};

void Circuit::add_gate(QuantumOperator const &gate){
    for (pos_t pos : gate.positions()){
        if (pos > qubit_num_) {
            throw "invalid position on quantum registers";
        }
    }
    gates_.push_back(gate);
}

 Circuit::Circuit(){};
 Circuit::Circuit(uint qubit_num)
 :
 qubit_num_(qubit_num), 
 gates_(0),
 max_targe_num_(0),
 measures_({})
 { }

 Circuit::Circuit(uint qubit_num, vector<QuantumOperator> const &gates, std::map<int, int> const& measures)
 :
 qubit_num_(qubit_num),
 gates_(gates),
 max_targe_num_(0),
 measures_(measures){}

Circuit::Circuit(py::object const&pycircuit, bool get_full_mat, bool reverse)
:
max_targe_num_(0)
{
    auto pygates = pycircuit.attr("operations");
    qubit_num_ = pycircuit.attr("num").cast<uint>();
    py::dict meas = pycircuit.attr("measures");
    // measures_ = py::list(meas.attr("keys")()).cast<vector<pos_t>>();
    for (auto item : meas){
        int qbit = item.first.cast<int>();
        int cbit = item.second.cast<int>();
        measures_[qbit] = cbit;
    }

    for (auto pygate_h : pygates){
        py::object pygate = py::reinterpret_borrow<py::object>(pygate_h);
        if (py::hasattr(pygate, "circuit")){
            auto wrap_circuit = Circuit(pygate.attr("circuit"), get_full_mat, reverse);
            for (auto op : wrap_circuit.gates()){
                gates_.push_back(op);
            }
        }
        else{
        QuantumOperator gate = from_pyops(pygate, get_full_mat, reverse);
        // check_operator(gate);
        if (gate){
            if (gate.targe_num() > max_targe_num_)
                max_targe_num_ = gate.targe_num();
            gates_.push_back(gate);
            }
        }        
    }
} 

