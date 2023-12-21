#pragma once

#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <Eigen/Core>
#include "types.hpp"
#include "util.h"

namespace py = pybind11;
using namespace pybind11::literals;

std::unordered_map<string, Opname> NAMED_GATE{Pair(creg), Pair(x), Pair(y), Pair(z), Pair(h), Pair(s), Pair(sdg), Pair(t), Pair(tdg), Pair(p), Pair(rx), Pair(ry), Pair(rz), Pair(cnot), Pair(cx), Pair(cz), Pair(cp), Pair(ccx), Pair(toffoli), Pair(rzz), Pair(measure)};

class QuantumOperator{
    protected:
        string name_;
        vector<pos_t> positions_;
        vector<double> paras_;
        uint control_num_;
        uint targe_num_;  
        bool diag_;
        bool real_;
        RowMatrixXcd targ_mat_;
        RowMatrixXcd full_mat_ ; 

    public:
        //Constructor
        QuantumOperator();
        QuantumOperator(string const name, vector<double> paras, vector<pos_t> const &control_qubits, vector<pos_t> const &targe_qubits, RowMatrixXcd const &mat = RowMatrixXcd(0, 0), RowMatrixXcd const &full_mat = RowMatrixXcd(0, 0), bool diag=false, bool real=false);

        QuantumOperator(string const name,vector<double> paras, vector<pos_t> const &positions, uint control_num, RowMatrixXcd const &mat=RowMatrixXcd(0, 0), RowMatrixXcd const &full_mat = RowMatrixXcd(0, 0), bool diag=false, bool real=false);

        QuantumOperator(string const name, vector<pos_t> const &positions, RowMatrixXcd const &mat=RowMatrixXcd(0, 0), RowMatrixXcd const &full_mat = RowMatrixXcd(0, 0));
        // virtual ~QuantumOperator();

        //data accessor
        string name() const {return name_;}
        vector<double> paras() const {return paras_;}
        bool has_control() const{return control_num_ == 0 ? false : true;}
        bool is_real() const{ return real_; }
        bool is_diag() const{ return diag_; }
        RowMatrixXcd targ_mat() const { return targ_mat_;}
        RowMatrixXcd full_mat() const {return full_mat_;}
        uint control_num() const { return control_num_; } 
        uint targe_num() const { return targe_num_; }
        vector<pos_t> positions() const { return positions_; }
        explicit operator bool() const {
            return !(name_ == "empty");
        }
};


QuantumOperator::QuantumOperator() : name_("empty"){ };

QuantumOperator::QuantumOperator(string const name, vector<pos_t> const &positions, RowMatrixXcd const &mat, RowMatrixXcd const &full_mat)
:
name_(name),
paras_(vector<double>(0)),
positions_(positions),
control_num_(-1),
targe_num_(0),
diag_(false),
real_(false),
targ_mat_(mat),
full_mat_(full_mat){ }

QuantumOperator::QuantumOperator(string const name, vector<double> paras, vector<pos_t> const &positions, uint control_num, RowMatrixXcd const &mat, RowMatrixXcd const &full_mat, bool diag, bool real)
:
name_(name),
paras_(paras),
positions_(positions),
control_num_(control_num),
targe_num_(positions.size()-control_num),
diag_(diag),
real_(real),
targ_mat_(mat), 
full_mat_(full_mat){ }

QuantumOperator::QuantumOperator(string const name, vector<double> paras, vector<pos_t> const &control_qubits, vector<pos_t> const &targe_qubits, RowMatrixXcd const &mat, RowMatrixXcd const &full_mat, bool diag, bool real)
:
name_(name),
paras_(paras),
diag_(diag),
real_(real),
targ_mat_(mat),
full_mat_(full_mat){
    positions_ = control_qubits;
    positions_.insert(positions_.end(), targe_qubits.begin(), targe_qubits.end());
    control_num_ = control_qubits.size();
    targe_num_ = targe_qubits.size();
}


class MeasuresOp: public QuantumOperator{
        private:
            std::vector<int> cbits_;
        public:
            MeasuresOp(){ };
            MeasuresOp(std::map<int, int> measures)
            {   name_ = "measure";
                for (auto it : measures){
                    positions_.push_back(it.first);
                    cbits_.push_back(it.second);
                }
            }
          
            std::vector<int> cbits(){ return cbits_;}
};

// Construct C++ operators from pygates
QuantumOperator from_pyops(py::object const &obj, bool get_full_mat=false, bool reverse=true){
    string name;
    vector<pos_t> positions;
    vector<double> paras;
    uint control_num = 0;
    RowMatrixXcd targ_mat;
    RowMatrixXcd full_mat;
    name = obj.attr("name").attr("lower")().cast<string>();
    if (!(name == "barrier" || name == "delay" || name == "id"))
    {   
        positions = obj.attr("pos").cast<vector<pos_t>>();
        paras = obj.attr("_paras").cast<vector<double>>();
        if (py::hasattr(obj, "ctrls")){
                control_num = py::len(obj.attr("ctrls"));
        }
    
        if (NAMED_GATE.count(name) == 0){
            if (py::hasattr(obj, "_targ_matrix")){
                targ_mat = obj.attr("_get_targ_matrix")("reverse_order"_a=true).cast<RowMatrixXcd>();
            }else{ //Single gate
                targ_mat = obj.attr("matrix").cast<RowMatrixXcd>();
            }
            
        }
        else{
            targ_mat = RowMatrixXcd(0, 0);
        }

        if (get_full_mat){
            full_mat = obj.attr("_get_raw_matrix")("reverse_order"_a=reverse).cast<RowMatrixXcd>();
        }else{
            full_mat = RowMatrixXcd(0, 0);
        }
        return QuantumOperator(name, paras, positions, control_num, std::move(targ_mat), std::move(full_mat));
    }
    
    else{
        return QuantumOperator();
    }
}

RowMatrixXcd extendMat(RowMatrixXcd const& targ_mat_, std::vector<pos_t> posv, int qnum){
    RowMatrixXcd mat = RowMatrixXcd::Zero(1<<qnum, 1<<qnum);
    auto posv_sorted = posv;
    sort(posv_sorted.begin(), posv_sorted.end());
    auto matsize = 1 << posv.size();
    std::vector<uint> targ_mask(matsize);
    for (size_t m = 0; m < matsize;m++){
        for (size_t j = 0; j < posv.size(); j++){
            if ((m>>j)&1){
                auto mask_pos = posv[j];
                targ_mask[m] |= 1ll<<mask_pos;
            }
        }
    }

    size_t rsize = 1<<(qnum-posv.size());
    for (size_t j = 0;j < rsize;j++){
        size_t i = j;
        for(size_t k=0;k < posv.size();k++){
            size_t _pos = posv_sorted[k];
            i = (i&((1ll<<_pos)-1)) | (i>>_pos<<_pos<<1);
        }
    
        for (size_t m1 = 0; m1 < matsize;m1++)
        for (size_t m2 = 0; m2 < matsize;m2++)
        {
            mat(i | targ_mask[m1], i | targ_mask[m2]) =
            targ_mat_(m1, m2);
        }
    }
    return mat;
}

QuantumOperator merge_operator(QuantumOperator const& op1, QuantumOperator const& op2){
    //TODO: use unitary matrix simulator 
    auto pos1 = op1.positions();
    auto pos2 = op2.positions();
    RowMatrixXcd extmat;
    if (pos1 == pos2){
        extmat = op2.full_mat() * op1.full_mat();
        return QuantumOperator("merged", {}, pos1, 0, extmat, extmat);
    }
    else{
        sort(pos1.begin(), pos1.end());
        sort(pos2.begin(), pos2.end());    
        if (std::includes(pos1.begin(), pos1.end(), pos2.begin(), pos2.end())){
            for (int j = 0;j < pos1.size();j++){
                for (int i2 = 0;i2 < pos2.size();++i2){
                    if (op1.positions()[j] == op2.positions()[i2]){
                        pos2[i2] = j;
                    }
                }
            }

            auto extmat2 = extendMat(op2.full_mat(), pos2, pos1.size());
            extmat = extmat2 * op1.full_mat();
            return QuantumOperator("merged", {}, op1.positions(), 0, extmat, extmat);
        }
        else if(std::includes(pos2.begin(), pos2.end(), pos1.begin(), pos1.end())){
             for (int j = 0;j < pos2.size();j++){
                for (int i1 = 0;i1 < pos1.size();++i1){
                    if (op2.positions()[j] == op1.positions()[i1]){
                        pos1[i1] = j;
                    }
                }
             }
            auto extmat1 = extendMat(op1.full_mat(), pos1, pos2.size());
            extmat = op2.full_mat() * extmat1;
            return QuantumOperator("merged", {}, op2.positions(), 0, extmat, extmat);
        }
        else{
            std::vector<pos_t> pos;
            set_union(pos1.begin(), pos1.end(), pos2.begin(), pos2.end(), back_inserter(pos));
           int qnum = pos.size();
            for (int j = 0;j < pos.size();j++){
                for (int i1 = 0;i1 < pos1.size();++i1){
                    if (pos[j] == op1.positions()[i1]){
                        pos1[i1] = j;
                    }
                }
                
                for (int i2 = 0;i2 < pos2.size();++i2){
                    if (pos[j] == op2.positions()[i2]){
                        pos2[i2] = j;
                    }
                }
            }
        
            auto extmat1 =  extendMat(op1.full_mat(), pos1, qnum);
            auto extmat2 =  extendMat(op2.full_mat(), pos2, qnum);
            extmat = extmat2 * extmat1;
            return QuantumOperator("merged", {}, pos, 0, extmat, extmat);
        }
    }
}
 

void check_operator(QuantumOperator &op){
    std::cout << "-------------" << std::endl;

    std::cout << "name: " << op.name() << std::endl;
    std::cout << "pos: ";
    Utils::printVector(op.positions());

    std::cout << "paras: ";
    Utils::printVector(op.paras());

    std::cout << "control number: ";
    std::cout << op.control_num() << std::endl;

    std::cout << "matrix: " << std::endl;
    std::cout << op.targ_mat() << std::endl;
    std::cout << "flatten matrix: " << std::endl;
    auto mat = op.targ_mat();
    // Eigen::Map<Eigen::RowVectorXcd> v1(mat.data(), mat.size());
    // std::cout << "v1: " << v1 << std::endl;
    auto matv = mat.data();
    for (auto i = 0;i < mat.size();i++){
        std::cout << matv[i] << " ";
    }
    std::cout << std::endl;

    std::cout << "full matrix: " << std::endl;
    std::cout << op.full_mat() << std::endl;
    std::cout << "-------------" << std::endl;
}

