#pragma once
#include <iostream>
#include <map>
#include <vector>
#include <queue>
#include <iqsimu/operators.hpp>
#include <iqsimu/circuit.hpp>
#include <set>
#include <algorithm>
#include <functional>


namespace DAG{

using Adjacents = std::map<int, std::map<int, int>>;

class DAGCircuit{
    private:
        int id_ = 0;
        int qubit_num_;
        std::map<int, QuantumOperator> nodes_;
        Adjacents adj;
        std::map<int, int> indegs;
        std::map<int, int> end_gates;
        std::map<int, int> begin_gates;
        MeasuresOp meas_;
    public:
        DAGCircuit(){ };
        DAGCircuit(Circuit const& circuit);
        int qubit_num(){ return qubit_num_; }
        void add_gate(QuantumOperator const &op){
            add_node(op);
            for (auto q : op.positions()){
                if (end_gates.count(q) > 0){
                    add_edge(end_gates[q], id_-1, q);
                }
                end_gates[q] = id_-1;
                if (begin_gates.count(q) == 0){
                    begin_gates[q] = id_-1;
                }
            }
        }

        void add_node(QuantumOperator const &op){
            nodes_[id_] = op;
            adj[id_] = std::map<int, int>();
            indegs[id_] = 0;
            id_++;
        }

        void add_edge(int src, int dest, int weight){
            if (adj[src].count(dest) == 0){
                adj[src][dest] = weight;
                indegs[dest] += 1;
            }
            else{
                adj[src][dest] = weight;
            }
        }

        void remove_gate(int ind)
        {   //TODO: handle begin and end
            auto behind = adj[ind];
            for (auto be : behind){
                indegs[be.first]--;
            }

            for (auto i : adj){
                auto front = i.second;
                for (auto fe : front){
                    if (fe.first == ind){
                        remove_edge(i.first, ind);
                        for (auto be : behind){
                            auto q = Utils::vectors_intersection(nodes_[i.first].positions(), nodes_[be.first].positions());
                            if (!q.empty()){
                                add_edge(i.first, be.first, q[0]);
                            }
                        }
                    }
                }
            } 

            nodes_.erase(ind);
            adj.erase(ind);
            indegs.erase(ind);
        }

        void remove_edge(int src, int dest){
            adj[src].erase(dest);
            indegs[dest]--;
        }

        bool are_adjacent(int n, int m){
            auto c = adj[n];
            auto q = c.find(m);
            if (q != c.end()){
                return true;
            }
            else{
                return false;
            }
        }

        void print(){
            for (auto i : nodes_){
                std::cout << i.first << "_" << i.second.name();
                if ( !adj[i.first].empty()){
                    std::cout << " {";
                    for (auto it : adj[i.first]){
                            std::cout<< " -" <<  it.second << "->" << it.first << "_" << nodes_[it.first].name();
                    }
                    std::cout << " }";
                }
                std::cout << " indeg: " << indegs[i.first] << std::endl;
                
            }
        }   

        std::map<int, QuantumOperator> nodes() const { return nodes_;}
        std::map<int, int> adjacents(int n) const { return adj.at(n);}
        int indeg(int n) const { return indegs.at(n); }
        int size() const { return nodes_.size(); }
        Adjacents adjacents() const { return adj;}
        void set_measures(MeasuresOp meas){
            meas_ = meas;
        }
        int add_measures(){
            add_gate(meas_);
            return id_-1;
        }
        Circuit to_circuit(){
            auto gates = topological_sort();
            std::map<int, int> measures;
            for (auto i = 0; i < meas_.positions().size();i++){
                measures[meas_.positions()[i]] = meas_.cbits()[i];
            }
            auto circuit = Circuit(qubit_num_, gates, measures);
            return circuit;
        }


        //algorithms
        std::vector<int> topological_order();
        std::vector<QuantumOperator>  topological_sort();
        void reverse();
        Adjacents DAGCircuit::reverse_adj();
        DAGCircuit find_connected(std::vector<int> svec);
        void DAGCircuit::absorb_gates();
        void merge_gates(int max_block_size);
        std::vector<int> merge_gate(int src, int dest, Adjacents & revadj);
       
};

DAGCircuit::DAGCircuit(Circuit const& circuit)
{   
    qubit_num_ = circuit.qubit_num();
    for (auto op : circuit.gates()){
        add_gate(op);
    }
    set_measures(MeasuresOp(circuit.measures()));
}

std::vector<int> DAGCircuit::topological_order(){
    std::queue<int> q;
    std::vector<int> order;
    auto indegs_ = indegs;
    for (auto it : indegs_){
        if (it.second == 0){
            q.push(it.first);
        }
    }

    while (!q.empty()){
        int u = q.front();
        // std::cout << "front" << u << std:: endl;
        q.pop();
        order.push_back(u);
        for (auto it : adj[u]){
            --indegs_[it.first];
            if (indegs_[it.first] == 0){
                q.push(it.first);
            }
        }
    }
    if (order.size() != indegs.size()){
        std::cout << "error: find cycle in dag" << std::endl;
        throw "error: find cycle in dag";
    }
    return order;
}

std::vector<QuantumOperator> DAGCircuit::topological_sort(){
    std::vector<QuantumOperator> oplist;
    auto order = topological_order();
    for (auto i : order){
        oplist.push_back(nodes_[i]);
    }
    return oplist;
}


void DAGCircuit::reverse(){
     auto adj_ = adj;
     for (auto it : indegs){
        indegs[it.first] = 0;
     }

     for (auto it: adj_){
        adj[it.first].clear();
     }

     for (auto i : adj_){
        for (auto it : i.second){
            add_edge(it.first, i.first, it.second);
        }
     }
     std::swap(begin_gates, end_gates);
};

Adjacents DAGCircuit::reverse_adj(){
    Adjacents revadj;
    for (auto i : adj){
        for (auto it : i.second){
            revadj[it.first][i.first] = it.second;
        }
    }
    return revadj;
}

DAGCircuit DAGCircuit::find_connected(std::vector<int> svec){
    DAGCircuit subDag;
    std::queue<int> q;
    std::map<int, int> visited;
    std::map<int, int> val_inds;

    for (auto i : svec){
        q.push(i);
    }

    int ind = 0;
    while (!q.empty()){
        int u = q.front();
        q.pop();
        if (visited.count(u) == 0){
            visited[u] = 1;
            subDag.add_node(nodes_[u]);
            val_inds[u] = ind;
            ++ind;
        }
        for (auto it : adj[u]){
            if (visited.count(it.first) == 0){
                visited[it.first] = 1;
                subDag.add_node(nodes_[it.first]);
                val_inds[it.first] = ind;
                ++ind;
                q.push(it.first);
            }
            subDag.add_edge(val_inds[u], val_inds[it.first], it.second);
        }
    }
    return subDag;
}

void DAGCircuit::absorb_gates(){
    std::queue<int> q;
    std::map<int, int> in_queue = indegs;
    for (auto it : in_queue){
        in_queue[it.first] = 0;
    }
    for (auto it : indegs){
        if (it.second == 0){
            q.push(it.first);
            in_queue[it.first]= 1;
        }
    }
    auto order = topological_order();
    std::map<int, int> order_dict;
    for (auto i = 0;i < order.size(); i++){
        order_dict[order[i]] = i;
    }

    while (!q.empty()){
        int u = q.front();
        // std::cout << "front" << u << std::endl;
        q.pop();
        in_queue[u] = 0;
        auto edges = adj[u];
        if (!edges.empty()){
        //sort the edges in order
            std::vector<int> sorted_edges;
            for (auto e : edges){
                sorted_edges.push_back(e.first);
            }
            sort(sorted_edges.begin(), sorted_edges.end(), [&order_dict](int x, int y){
                return order_dict[x] < order_dict[y];
                });
            std::vector<uint> excluded_pos;
            for (auto next : sorted_edges){  
                //merge gate and add gate to the queue
                auto op1 = nodes_[u];
                auto op2 = nodes_[next];
                auto can_merge = Utils::vectors_intersection(op2.positions(), excluded_pos).empty(); 
                if (Utils::vectors_includes(op1.positions(), op2.positions())){
                    if (can_merge){
                        // std::cout << "merge: " << u << " " << next << std::endl;
                        nodes_[u] = merge_operator(op1, op2);
                        if (in_queue[u] == 0){
                            q.push(u);
                            in_queue[u] = 1;
                        }
                        remove_gate(next);
                    }
                }
                else if (Utils::vectors_includes(op2.positions(), op1.positions())){
                    if (can_merge){
                        // std::cout << "merge: " << next << " " << u << std::endl;
                        nodes_[next] = merge_operator(op1, op2);
                        if (in_queue[next] == 0){
                            q.push(next);
                            in_queue[next] = 1;
                        }
                        remove_gate(u);
                    }
                }
                else{
                    // std::cout << "not include : " << next << " " << u << std::endl;
                    if (in_queue[next] == 0){
                        q.push(next);
                        in_queue[next] = 1;
                    }
                    for (uint pos : op2.positions()){
                        excluded_pos.push_back(pos);
                    }
                }
            }
        } 
    }
}

std::vector<int>  DAGCircuit::merge_gate(int u, int v, Adjacents & revadj){
    auto op1 = nodes_[u];
    auto op2 = nodes_[v];
    
    nodes_[u] = merge_operator(op1, op2);
  
    //update behind
    std::vector<int> newedges;
    for (auto e : adj[v]){
        indegs[e.first]--;
        revadj[e.first].erase(v);
        if (e.first != u){
            add_edge(u, e.first, e.second);
            revadj[e.first][u] = e.second;
            newedges.push_back(e.first);
            // std::cout << "add_behind " << u <<  "->" << e.first << std::endl;
        }
    }

    //update front
    for (auto e : revadj[v]){
        remove_edge(e.first, v);
        if (e.first != u){
            add_edge(e.first, u, e.second);
            revadj[u][e.first] = e.second;
            // std::cout << "add_front " << e.first << "->" << u << std::endl;
        }
    }


    nodes_.erase(v);
    adj.erase(v);
    indegs.erase(v);
    revadj.erase(v);

    return newedges;
}   

bool if_obstacled(int src, int dest, std::map<int, int> & order_dict, Adjacents & revadj, std::set<int> excluded_nodes){
        bool obstacled = false;
        std::queue<int> q;
        std::map<int, int> visited;
        q.push(src);
        while (!q.empty()){
            int u = q.front();
            q.pop();
            if (visited.count(u) == 0){
                visited[u] = 1;
                if (excluded_nodes.count(u)){
                    obstacled = true;
                    break;
                }
            }
            for (auto it : revadj[u]){
                if (visited.count(it.first) == 0){
                    visited[it.first] = 1;
                    if (excluded_nodes.count(it.first)){
                        obstacled = true;
                        break;
                    }
                    if (order_dict[it.first] > order_dict[dest]){
                       q.push(it.first);
                    }                                  
                }
            }
        }
        return obstacled;
}


void DAGCircuit::merge_gates(int max_block_size){
    std::queue<int> q;
    std::map<int, int> in_queue = indegs;
    for (auto it : in_queue){
        in_queue[it.first] = 0;
    }
    for (auto it : indegs){
        if (it.second == 0){
            q.push(it.first);
            in_queue[it.first]= 1;
        }
    }
    
    auto revadj = reverse_adj();
    while (!q.empty()){
        int u = q.front();
        // std::cout << "front " << u << std::endl;
        q.pop();
        in_queue[u] = 0;
        auto edges = adj[u];
        if (!edges.empty()){
            //sort the edges in order
            auto order = topological_order();
            std::map<int, int> order_dict;
            for (auto i = 0;i < order.size(); i++){
                order_dict[order[i]] = i;
            }
            auto topoComp = [&order_dict](int x, int y){
                return order_dict[x] > order_dict[y];
            };
            std::vector<int> sorted_edges;
            for (auto e : edges){
                sorted_edges.push_back(e.first);
            }
            std::make_heap(sorted_edges.begin(), sorted_edges.end(), topoComp);
            std::set<int> excluded_nodes;
            bool excluded = false;

            while(!sorted_edges.empty()){
                // std::cout << "edges" << std::endl;
                // Utils::printVector(sorted_edges);
                auto next = sorted_edges.front();
                // std::cout << "-------------edge: " << u << "->" << next << "--------" << std::endl;
                auto op1 = nodes_[u];
                auto op2 = nodes_[next];
                auto op2_pos = op2.positions();
                // Utils::printVector(op1.positions());
                // Utils::printVector(op2_pos);
                
                bool exceed = Utils::vectors_union(op1.positions(), op2_pos).size() > max_block_size;
                bool obstacle = false;
                
                if (excluded){
                    //find obstacle
                    obstacle = if_obstacled(next, u, order_dict, revadj, excluded_nodes);
                }

                if (!exceed && !obstacle){
                    //merge gate and add gate to the queue
                    // std::cout << "merge: " << u << "->" << next << std::endl;
                    auto newedges = merge_gate(u, next, revadj);
                    std::pop_heap(sorted_edges.begin(), sorted_edges.end(), topoComp);
                    sorted_edges.pop_back();
                    
                    for (auto ne : newedges){
                        for (auto iter = sorted_edges.begin();iter !=sorted_edges.end();){
                            if (*iter == ne){
                                iter = sorted_edges.erase(iter);
                            }
                            else{ iter++ ;}
                        }
                        sorted_edges.push_back(ne);
                        std::push_heap(sorted_edges.begin(), sorted_edges.end(), topoComp);
                    }
                }
                else{
                    // std::cout << "exclude " << next << std::endl;
                    if (in_queue[next] == 0){
                        q.push(next);
                        in_queue[next] = 1;
                    }
                    excluded_nodes.insert(next);
                    excluded = true;
                    std::pop_heap(sorted_edges.begin(), sorted_edges.end(), topoComp);
                    sorted_edges.pop_back();
                }
            }
        }
    }
}

 //TODO:merge unconnected small blocks

} //name space DAG