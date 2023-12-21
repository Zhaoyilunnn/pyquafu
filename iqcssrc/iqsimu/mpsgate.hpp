#pragma once
#include "types.hpp"
#include "itensor/all.h"

using namespace itensor;

namespace IQMPS{

class iTenGate{
    protected:
        Index _q1;
        Index _q2;
        ITensor _gateT;
        int _qbit_num;
    public:
        iTenGate(){ };
        iTenGate(Index q1, RowMatrixXcd const & mat);
        iTenGate(Index q1, Index q2, RowMatrixXcd const &mat);
        virtual ~iTenGate(){ };

        Index const q1(){ return _q1; }
        Index const q2(){ return _q2; }
        int qubitNum(){ return _qbit_num;}
        ITensor gateT(){ return _gateT; }
        void printGate();
};

iTenGate::iTenGate(Index q1, RowMatrixXcd const &mat) :
_q1(q1),
_qbit_num(1)
{
    _gateT = ITensor(_q1, prime(q1));
    #pragma unroll
    for (int i=0;i<2;i++)
    for (int j=0;j<2;j++){
        _gateT.set(_q1=i+1, prime(_q1)=j+1, mat(j, i));
    }
}

iTenGate::iTenGate(Index q1, Index q2, RowMatrixXcd const &mat) :
_q1(q1),
_q2(q2),
_qbit_num(2)
{
    _gateT = ITensor(_q1, prime(q1), _q2, prime(_q2));
    #pragma unroll
    for (int i=0;i<2;i++)
    for (int j=0;j<2;j++)
    for (int pi=0;pi<2;pi++)
    for (int pj=0;pj<2;pj++){
        _gateT.set(_q1=i+1, prime(_q1)=pi+1, _q2=j+1, prime(_q2)=pj+1, mat(pi*2+pj, i*2+j));
    }
}


void iTenGate::printGate()
{   
    PrintData(_gateT);
}

class MPSGate : public iTenGate{
    private:
        int _pos1;
        int _pos2;

    public:
        MPSGate(){ };
        MPSGate(SiteSet &site, int pos1, RowMatrixXcd const &mat);
        MPSGate(SiteSet &site, int pos1, int pos2, RowMatrixXcd const &mat);
        virtual ~MPSGate(){ };

        int const pos1(){ return _pos1; }
        int const pos2(){ return _pos2; }
};

MPSGate::MPSGate(SiteSet &site, int pos1, RowMatrixXcd const &mat)
:
_pos1(pos1),
_pos2(-1),
iTenGate(site(pos1), mat)
{  }

MPSGate::MPSGate(SiteSet &site, int pos1, int pos2, RowMatrixXcd const &mat)
:
_pos1(pos1),
_pos2(pos2),
iTenGate(site(pos1), site(pos2), mat)
{   if (pos1 == pos2){
        throw "invalid bits";
    } 
}

double apply_mpsgate(MPSGate & g, MPS & psi, Args args = Args::global()){
    double truncerr = 0.;
    if (g.qubitNum() == 1) //single qubit gate
    {   
        int i = g.pos1();
        auto A = psi(i) * g.gateT();
        A.noPrime();
        psi.set(i, A);
    }
    else if (g.qubitNum() == 2){
        int i1, i2;
        if (g.pos1() > g.pos2()){
            i1 = g.pos2();
            i2 = g.pos1();
        }else{
            i1 = g.pos1();
            i2 = g.pos2();
        }
        psi.position(i1);
        if (i2 > i1 + 1){
            for (auto j = i2; j-1>=i1+1;--j)
                {
                    auto temp = psi(j-1)*psi(j);
                    auto inds = IndexSet(uniqueIndex(psi(j-1), psi(j), "Link"), findIndex(psi(j), "Site"));
                    auto [U, S, V] = svd(temp, inds, args);
                    psi.set(j-1, U*S);
                    psi.set(j, V);
                    psi.ref(j-1) /= norm(psi(j-1));
                }

            auto AA = psi(i1)*psi(i1+1)*g.gateT();
            AA.noPrime();
            ITensor d;
            svd(AA, psi.ref(i1), d, psi.ref(i1+1), args);
            truncerr = std::max(sqr(norm(psi.ref(i1)*d*psi.ref(i1+1) - AA)/norm(AA)), truncerr);
            psi.ref(i1+1) *= d;
            psi.ref(i1+1) /= norm(psi(i1+1));
            //swap back
            for (auto j=i1+1; j+1<=i2;++j)
            {
                auto temp = psi(j)*psi(j+1);
                auto  inds = IndexSet(uniqueIndex(psi(j), psi(j+1), "Link"), findIndex(psi(j+1), "Site"));
                auto [U, S, V] = svd(temp, inds, args);
                psi.set(j, U);
                psi.set(j+1, S*V);
                psi.ref(j+1) /= norm(psi(j+1));
            }
        }
        else{
            auto AA = psi(i1)*psi(i1+1)*g.gateT();
            AA.noPrime();
            ITensor d;
            svd(AA, psi.ref(i1), d, psi.ref(i1+1), args);
            truncerr = std::max(sqr(norm(psi.ref(i1)*d*psi.ref(i1+1) - AA)/norm(AA)), truncerr);
            psi.ref(i1+1) *= d;
            psi.ref(i1+1) /= norm(psi(i1+1));
        }
    }
    else{
        throw "invalid bit";
    }
    return truncerr;
}

}// namespace IQMPS


