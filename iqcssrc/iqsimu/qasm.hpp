#pragma once

#include <unordered_map>
#include "types.hpp"
#include "util.h"

using Utils::split_string;
using Utils::find_numbers;

typedef struct{
    string name;
    vector<pos_t> positions;
    vector<double> params; 

    void print_info(){
        std::cout << "name " << name << std::endl;
        std::cout << "positions: ";
        for (auto pos : positions){
            std::cout << pos << " ";
        }
        std::cout << std::endl;

        if (params.size() > 0){
        printf("parameters: ");
            for (auto para : params){
                printf("%.6f ", para);
            }
        }
        printf("\n");
        printf("-----\n");
    }

} Operation;

Operation compile_line(string const& line){
    auto operation_qbits = split_string(line, ' ', 1);
    auto operation = operation_qbits[0];
    auto qbits = operation_qbits[1];
    auto positions = find_numbers<pos_t>(qbits);
    auto opname_params = split_string(operation, '(', 1);
    auto opname = opname_params[0]; 
    vector<double> params;
    if (opname_params.size() > 1){
        params = find_numbers<double>(opname_params[1]);
    }

    return Operation{opname, positions, params};
}