#pragma once

#include "circuit.hpp"
#include "mpsgate.hpp"
#include "statevector.hpp"
#include "optimize/reduce.hpp"
#include <string>
#include <unordered_map>
using ITallocator = itensor::uninitialized_allocator<std::complex<double>>;

auto extractCplx = [](Dense<Cplx> const& d){
    return d.store;
};

namespace IQMPS{

struct MpsResult{
    std::unordered_map<std::string, int> counts;
    std::vector<Cplx, ITallocator> ITVec;
    double max_truncerr = 0;
    std::vector<double> paulis_expects;
};

ITensor getStateVec(MPS &psi){
    ITensor Wf = psi(1);
    for (int i=2;i <= psi.length();i++){
        Wf = Wf * psi(i);
    }
    return Wf;
}

//measureMPS
std::unordered_map<std::string, int> measureMPS(MPS & psi, std::map<int, int> meas, int shots){
    std::unordered_map<std::string, int> counts;
    auto sites = siteInds(psi);
    auto N = length(psi);
    auto m_start = (meas.begin()->first)+1;
    auto end = meas.end();
    end--;
    auto m_end = (end->first)+1;
    psi.position(m_start);

    // OMP Parallel need OPENBLAS compiled with openmp
    int num_procs = omp_get_num_procs();
    int threads = 2*num_procs-2;
    if (threads > shots){threads = shots;}
#pragma omp parallel for num_threads(threads)
    for (int s = 0 ;s < shots;s++){
        std::string bitstring(meas.size(), '0');
        Index sj = sites(m_start);
        auto PUp = ITensor(sj, prime(sj));
        auto PDown = ITensor(sj,prime(sj));
        PUp.set(1, 1, 1.0);
        PDown.set(2, 2, 1.0);
        Real prob_last = 1.0;
        auto rl = commonIndex(psi(m_start), psi(m_start+1));
        auto E = psi(m_start) * PUp * dag(prime(prime(psi(m_start), "Site"), rl));
        auto Eup = E;
        Real prob_up = eltC(E * delta(rl, prime(rl))).real();
        prob_last = prob_up;

        if(Global::random() > prob_up){
            E = psi(m_start) * PDown * dag(prime(prime(psi(m_start), "Site"), rl));
            prob_last = 1 - prob_up;
            bitstring.replace(meas.begin()->second, 1, "1");
        }
        for (int m = m_start + 1;m <= m_end; m++){
            if (meas.count(m-1) == 0)//no measure
            {   
                E = E * psi(m) * dag(prime(psi(m), linkInds(psi, m)));
                if (m < N ){
                    rl = commonIndex(psi(m), psi(m+1));
                }
            }
            else{//perform measure
                Index sj = sites(m);
                auto PUp = ITensor(sj, prime(sj));
                auto PDown = ITensor(sj,prime(sj));
                PUp.set(1, 1, 1.0);
                PDown.set(2, 2, 1.0);
                Eup = E * psi(m)* PUp * dag(prime(psi(m)));
                if (m < N ){
                    rl = commonIndex(psi(m), psi(m+1));
                    prob_up = eltC(Eup * delta(rl, prime(rl))).real();
                }else{
                    prob_up = eltC(Eup).real();
                }
            
                if(Global::random() > prob_up / prob_last){
                    E =  E * psi(m) * PDown * dag(prime(psi(m)));
                    prob_last = prob_last -  prob_up;
                    bitstring.replace(meas[m-1], 1, "1");
                }
                else{
                    E = Eup;
                    prob_last = prob_up;
                }
            } 
        }
        if (counts.count(bitstring) >  0){
#pragma omp atomic
            counts[bitstring] += 1;
        }else{
#pragma omp critical
            counts[bitstring] = 1;
        }
    }
    return counts;
}


ITensor pauliTensor(Index const& sj, char pauli){
    auto localop = ITensor(sj, prime(sj));
    switch(pauli){
        case 'X':{
            localop.set(1, 2, 1.0);
            localop.set(2, 1, 1.0);
            break;
        }
        case 'Y': {
            localop.set(1, 2, imag_I);
            localop.set(2, 1, -imag_I);
            break;
        }
        case 'Z': {
            localop.set(1, 1, 1.0);
            localop.set(2, 2, -1.0);
            break;
        }
    }
    return localop;
}

double expect_pauli(MPS & psi, std::map<int, char> & pauli){
    auto sites = siteInds(psi);
    auto N = length(psi);
    auto m_start = (pauli.begin()->first)+1;
    auto end = pauli.end();
    end--;
    auto m_end = (end->first)+1;
    psi.position(m_start);
    Index sj0 = sites(m_start);
    auto op0 =  pauliTensor(sj0, pauli[m_start-1]);
    Index rl;
    ITensor E;
    if (m_start < N){
        rl = commonIndex(psi(m_start), psi(m_start+1));
        E = psi(m_start) * op0 * dag(prime(prime(psi(m_start), "Site"), rl));
    }
    else{
         return  eltC(psi(m_start) * op0 * dag(prime(psi(m_start), "Site"))).real();
    }
    for (int m = m_start + 1;m <= m_end;m++){
        if (pauli.count(m-1) == 0){
             E = E * psi(m) * dag(prime(psi(m), linkInds(psi, m)));
             if (m < N){
                rl = commonIndex(psi(m), psi(m+1));
             }
        }
        else{
            Index sj = sites(m);
            auto op =  pauliTensor(sj, pauli[m-1]);
            E = E * psi(m) * op * dag(prime(psi(m)));
            if (m < N){
                rl = commonIndex(psi(m), psi(m+1));
            }
            else{
                return eltC(E).real();
            }
        }
    }
    return eltC(E * delta(rl, prime(rl))).real();
}

std::vector<double> expect_paulis(MPS & psi,std::vector<std::map<int, char>> & paulis_maps){
    std::vector<double> expects;
    for (auto pauli_map : paulis_maps){
        auto res = expect_pauli(psi, pauli_map);
        expects.push_back(res);
    }
    return expects;
}

MpsResult simulate_mps(Circuit & circuit, std::vector<std::map<int, char>> & paulis_maps, double cutoff=1e-16, int shots=0,  bool save_statevec=false)
{
    MpsResult res;
    auto N = circuit.qubit_num();
    auto meas = circuit.measures(); 
    bool no_obs = paulis_maps.empty();
    //reduce
    if (no_obs && (!save_statevec || N > 26)){
        circuit = reduce_by_measure(circuit);
    }
    //initialize mps
    auto site = SpinHalf(N, {"ConserveQNs=",false});
    MPS psi = MPS(site);
    auto l = linkInds(psi);
    psi.ref(1) = setElt(l(1)=1, site(1)=1);
    for (auto i : range(2, N)){
        psi.ref(i) = setElt(site(i)(1), l(i-1)(1), l(i)(1));
    }
    psi.ref(N) = setElt(l(N-1)(1), site(N)(1));

    //apply circuit
    for (auto op : circuit.gates())
    {   
        MPSGate opT;
        if (op.positions().size() == 1){
            opT = MPSGate(site, op.positions()[0]+1, op.full_mat());
        }else if(op.positions().size() == 2){
            opT = MPSGate(site, op.positions()[0]+1, op.positions()[1]+1, op.full_mat()); 
        }
        else{
            throw "not support for gate act on more than 2 qubits";
        }
        auto truncerr = apply_mpsgate(opT, psi, {"Cutoff=", cutoff});  
        res.max_truncerr = std::max(res.max_truncerr, truncerr);      
    }
    if (shots > 0){
        auto counts = measureMPS(psi, meas, shots);
        res.counts = counts;
    }
    
    if (!no_obs){
        res.paulis_expects = expect_paulis(psi, paulis_maps);
    }

    if (save_statevec){
        if (N <= 26){
            ITensor Wf = getStateVec(psi);
            Wf = Wf * Cplx(1, 1E-25);
            auto wf = applyFunc(extractCplx, Wf.store());
            res.ITVec = std::move(wf);
        }
        else{
            std::cout << "warning: qubit number larger than 26, no statevector saved \n" << std::endl;
        }
    }
    return res;
}

}// namespace IQMPS